
#library(tidyverse)
library(purrr)
library(dplyr)
library(readr)
library(magrittr)
library(MASS)
library(seqinr)


#setwd("/Users/ahrnee-adm/dev/R/workspace/SpectroX/tests/testthat/")
#source("./R/SpectroX.R")

TESTFILE = "../../inst/testData/msms_test.txt"
TARGETSFILE = "../../inst/testData/proteins_acs.txt"
TMPXLS = paste0(tempdir(),"/tmp.xls")
TMPPDF = paste0(tempdir(),"/tmp.pdf")

TB = parseMaxQuantMSMS(TESTFILE)
IRTMODEL = getIRTModel(subset(TB, isIRT))
UPFASTAFILE = "../../inst/testData/protein_db1.fasta"
### read protein db
UPPROTEINDB <- read.fasta(UPFASTAFILE,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
PEPTIDEDF = digestProteome(UPPROTEINDB, peptideLengthRange = c(6,21), dispProgressBar = F)
SPECTRALLIBRARY = createSpectralLibrary(TB
                                        , minFragNb= 3
                                        , minNbTransitions= 4
                                        , minBasePeakFraction = 0.1
                                        , includeNL = TRUE
)

# unit tests
testParseMaxQuantMSMS = function(){

  cat("--- testParseMaxQuantMSMS: --- \n")

  stopifnot( !(parseMaxQuantMSMS(TESTFILE, label = "light")$Modifications %>% grepl("Arg10|Lys8",.)) )
  stopifnot( sum(parseMaxQuantMSMS(TESTFILE, label = "heavy", keepBestSpectrumOnly=F, selectedPTMRegExp = "Oxid")$Modifications %>% grepl("Arg10|Lys8",.)) == 39 )
  tb = parseMaxQuantMSMS(TESTFILE, pepCutoff = 0.01, filterNonExclusivePeptides = T,chargeState = 3, selectedPTMRegExp = "Oxid")
  stopifnot(tb$PEP < 0.01)
  stopifnot(!grepl(";",tb$Proteins))
  stopifnot(tb$Charge ==  3 |  grepl("Biognosys",tb$Proteins))
  stopifnot( sum( grepl("Oxid",tb$Modifications) ) > 0  )
  stopifnot( sum(rownames(IRTPEPTIDES) %in% tb$Sequence) == 9  )
  # no oxidation
  tb2 = parseMaxQuantMSMS(TESTFILE, pepCutoff = 0.01, filterNonExclusivePeptides = F)
  stopifnot( sum( grepl("Oxid",tb2$Modifications) ) == 0  )
  # do not filter non exclusive
  stopifnot(sum(grepl(";",tb2$Proteins)) > 0)
  # no contaminants
  stopifnot( !grepl("CON_",tb2$Proteins) | grepl("Biognosys",tb2$Proteins) )

  # target peptides
  tbTargetPep =  parseMaxQuantMSMS(TESTFILE, pepCutoff = 0.01, filterNonExclusivePeptides = T, targetPeptides = c("AAVPSGASTGIYEALELR","AFMNNK","AATFGLILDDVSLTHLTFGK") )
  stopifnot(nrow(subset(tbTargetPep, !isIRT)) == 2)

  # target proteins
  tbTargetProteins =  parseMaxQuantMSMS(TESTFILE, pepCutoff = 0.01, filterNonExclusivePeptides = F, targetProteins = c("P67778","P62259","Q9D0K2"), maxMissedCleavages = 2, keepBestSpectrumOnly=F )
  stopifnot(nrow(tbTargetProteins) == 24)

  # keepBestSpectrumOnly
  stopifnot(table(paste(tb$Sequence,tb$Modifications)) %>% max ==1)
  stopifnot((parseMaxQuantMSMS(TESTFILE, pepCutoff = 1, maxMissedCleavages = 2, keepBestSpectrumOnly = F ) %>% nrow) > nrow(tb))

  # mc
  tbMC = parseMaxQuantMSMS(TESTFILE, maxMissedCleavages = 2)
  stopifnot( tbMC$`Missed cleavages` %>% max == 2)

  # ignoreArgLysIsoLabel = F
  stopifnot( parseMaxQuantMSMS(TESTFILE, pepCutoff = 1, maxMissedCleavages = 2, ignoreArgLysIsoLabel=F )$Modifications %>% unique() == "Unmodified")

  # pep length
  tbPepLen = parseMaxQuantMSMS(TESTFILE,pepLength = c(10,Inf))
  stopifnot(nchar(tbPepLen$Sequence)> 10 |  grepl("Biognosys",tbPepLen$Proteins) )
  cat("--- testParseMaxQuantMSMS: PASS ALL TEST --- \n")

}

testGetSearchedModifications = function(){
  cat("--- testGetSearchedModifications: --- \n")

  #stopifnot( c("Arg10","Lys8","Oxidation (M)") %in%  getSearchedModifications(TB))
  stopifnot( c("Oxidation (M)") %in%  getSearchedModifications(TB))
  cat("--- testGetSearchedModifications: PASS ALL TEST --- \n")
}

testCreateLibrarySpectrum <- function(){

  cat("--- testCreateLibrarySpectrum: --- \n")
  annotSpec = createLibrarySpectrum(TB[7,])
  stopifnot(dim(annotSpec) == c(17,23) )
  stopifnot( c("a","b","y") == levels(annotSpec$ionType) )
  stopifnot(annotSpec$adjustedIntensity > annotSpec$intensity    )
  stopifnot(annotSpec$precCharge == TB[7,]$Charge)
  stopifnot(annotSpec$precApexIntensity >=    annotSpec$precIntensity)

  cat("--- testCreateLibrarySpectrum: PASS ALL TEST --- \n")
}


tesGetIRTModel <- function(){

  cat("--- tesGetIRTModel: --- \n")
  stopifnot(dim(IRTMODEL$irtPeptides) == c(10,3)  )
  stopifnot(coef(IRTMODEL$fit) %>% round == c(-89,3)  )
  cat("--- tesGetIRTModel: PASS ALL TEST --- \n")
}

testGetEmpiricalIRT <- function(){

  cat("--- testGetEmpiricalIRT: --- \n")
  tb = TB
  tb$iRT =    getEmpiricalIRT(tb = TB,fit = IRTMODEL$fit )

  #stopifnot( mean(  tb$iRT - IRTPEPTIDES[match(tb$Sequence, rownames(IRTPEPTIDES)),], na.rm=T ) %>% round == 0)
  stopifnot(!is.na(createLibrarySpectrum(tb[9,])$iRT) )
  cat("--- testGetEmpiricalIRT: PASS ALL TEST --- \n")
}

testGetFragmentSequence <- function(){
  cat("--- testGetFragmentSequence: --- \n")
  stopifnot(getFragmentSequence(peptide="PEPTIDEK",ionType="b",fragmentNb=3) == "PEP")
  stopifnot(getFragmentSequence(peptide="PEPTIDEK",ionType="a",fragmentNb=1) == "P")
  stopifnot(getFragmentSequence(peptide="PEPTIDEK",ionType="y",fragmentNb=3) == "DEK")
  stopifnot(getFragmentSequence(peptide="PEPTIDEK",ionType="x",fragmentNb=1) == "K")
  cat("--- testGetFragmentSequence:  PASS ALL TEST --- \n")
}

# testGetLabelMzShift <- function(){
#
#   cat("--- testGetLabelMzShift: --- \n")
#
#   ### ONLY LABELLED AT THE TERMINI
#
#   stopifnot(getLabelMzShift(aaSeq = "PETIDEK", charge=1, isHeavy=T) == -8.014199)
#   stopifnot(getLabelMzShift(aaSeq = "PEPTIDER", charge=1, isHeavy=T) == -10.008269)
#   stopifnot(getLabelMzShift(aaSeq = "PEKPTIDEK", charge=1, isHeavy=F) == 8.014199)
#   stopifnot(getLabelMzShift(aaSeq = "PERPTIDER", charge=2, isHeavy=F) == 10.008269/2)
#
#   cat("--- testGetLabelMzShift:  PASS ALL TEST --- \n")
#
# }

# testGetComplementaryIsotopeAnnotatedSpectrum <- function(){
#
#   cat("--- testGetComplementaryIsotopeAnnotatedSpectrum: --- \n")
#
#   annotSpec = createLibrarySpectrum(TB[9,])
#   compAnnotSpec = getComplementaryLabelSpectrum(annotSpec)
#
#   stopifnot(subset(annotSpec, ionType == "y")$mz < subset(compAnnotSpec, ionType == "y")$mz)
#   stopifnot(subset(annotSpec, ionType == "b")$mz == subset(compAnnotSpec, ionType == "b")$mz)
#   stopifnot(subset(annotSpec, ionType == "a")$mz == subset(compAnnotSpec, ionType == "a")$mz)
#   stopifnot(compAnnotSpec$isHeavy == !annotSpec$isHeavy)
#
#   cat("--- testGetComplementaryIsotopeAnnotatedSpectrum:  PASS ALL TEST --- \n")
#
# }

testGetPeptides <- function(){

  cat("--- testGetPeptides: --- \n")

  proteinSeq1 <- "MSAGSSCSQTPSRAIPTRRVALGDGVQLPPGDYSTTPGGTLFSTTPGGTRIIYDRKFLMECRNSPVAKTPPKDLPAIPGVTSPTSDEPPMQASQSQLPSSPEDKRAGGEESQFEMDI"
  proteinSeq2 <- "MVKKSRRRGAAQWAAVRAQAGLTATDENEDDLGLPPSPGDSSYYQDQVDEFHEARSRAVLAKGWNEVESGEEDGDEEEE"
  proteinSeq3 <- "MSERMRPVVVDLPTSASSSMKVNG"
  proteinSeq4 <- "RKR"

  stopifnot(length(getPeptides(proteinSeq3)) == 3)
  stopifnot(paste(getPeptides(proteinSeq3),collapse="") == proteinSeq3)
  stopifnot(length(getPeptides(proteinSeq3,proteaseRegExp="K(?!P)")) == 2)

  stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=0)) == 3)
  stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=1)) == 5)
  stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=2)) == 6)

  cat("--- testGetPeptides:  PASS ALL TEST --- \n")

}

testDigestProteome = function(){

  cat("--- testDigestProteome: --- \n")

  stopifnot( nrow(PEPTIDEDF) == 27)
  # get exclusive peptides only i.e. peptides mapping to a single protein
  stopifnot((PEPTIDEDF %>% nrow) == (table(PEPTIDEDF$peptide)[table(PEPTIDEDF$peptide) == 1 ] %>% length))
  stopifnot( nrow(digestProteome(UPPROTEINDB, peptideLengthRange = c(6,21), dispProgressBar = F, nbMiscleavages = 2 )) > 3*nrow(PEPTIDEDF) )
  stopifnot(nchar(digestProteome(UPPROTEINDB, peptideLengthRange = c(3,3), dispProgressBar = F)$peptide %>% as.character ) == 3)
  stopifnot(digestProteome(UPPROTEINDB, exclusivePeptides = T, dispProgressBar = F, trimAC=T )$protein %in% c("P62258","Q9JII5"))
  stopifnot(digestProteome(UPPROTEINDB, exclusivePeptides = F, dispProgressBar = F, trimAC=T )$protein %in% c("P62258", "P31946","P31946-2","Q9JII5"))
  stopifnot(digestProteome(UPPROTEINDB, exclusivePeptides = F, dispProgressBar = F, trimAC=F )$protein %>% grepl("sp",.))
  cat("--- testDigestProteome:  PASS ALL TEST --- \n")

}

testProteotypicPeptideExport = function(){

  cat("--- testProteotypicPeptideExport: --- \n")

  specLib = createLibrarySpectrum(TB[9,])
  specLib$rankingMetric = 1
  specLib$adjustedIntensitySum = 1
  specLib$psmScore = 1
  specLib$precApexIntensity = 1
  proteotypicPeptideExport(spectralLibrary = specLib
                           ,targetProteins = c(TB$Proteins, PEPTIDEDF$protein %>% unique())
                           , nbPeptidesPerProtein = 6
                           , theoPeptides =PEPTIDEDF
                           , outFile= paste0(tempdir(),"/bla.xls")
  )

  cat("--- testProteotypicPeptideExport:  PASS ALL TEST --- \n")

}

testCreateSpectralLibrary = function(){

  cat("--- testCreateSpectralLibrary: --- \n")

  arrange(subset(TB,Sequence == "GTFIIDPAAVIR") %>% createLibrarySpectrum(), intensity)[1:5,]
 # AFSEGQITR
  arrange(subset(TB,Sequence == "AFSEGQITR") %>% createLibrarySpectrum(), intensity)[1:5,]

  stopifnot(SPECTRALLIBRARY$fragmentNb >= 3)
  stopifnot((subset(SPECTRALLIBRARY,SPECTRALLIBRARY$mqResIdx == unique(SPECTRALLIBRARY$mqResIdx)[1])$precMz %>% unique %>% length) == 1)
  stopifnot(SPECTRALLIBRARY$rankingMetric == SPECTRALLIBRARY$adjustedIntensitySum)

  spectralLibrary2 = createSpectralLibrary(TB
                                           , minFragNb= 3
                                           , minNbTransitions= 3
                                           , maxNbTransitions = 5
                                           , minBasePeakFraction = 0.1
                                           , ionTypeFilter = c("b")
                                           , includeNL = F
                                           , rankingMetric = "psmScore"
  )

  # rankingMetric
  stopifnot(spectralLibrary2$rankingMetric == spectralLibrary2$psmScore)
  # minNbTransitions and maxNbTransitions
  transitionsPerSpec =spectralLibrary2  %>%  split(.$mqResIdx) %>%
    map(~ nrow(.x)) %>% unlist
  stopifnot(min(transitionsPerSpec) == 3 & max(transitionsPerSpec) == 5)

  stopifnot(spectralLibrary2$ionType == "b")
  cat("--- testCreateSpectralLibrary:  PASS ALL TEST --- \n")

}


testCreateComplementaryIsotopeLibrary = function(){

  cat("--- testCreateComplementaryIsotopeLibrary:  --- \n")

  compSL = createComplementaryIsotopeLibrary(SPECTRALLIBRARY)

  # light -> heavy
  stopifnot(sum(!compSL$isHeavy)  == sum(SPECTRALLIBRARY$isHeavy))
  # all precMz should be different
  stopifnot(compSL$precMz != SPECTRALLIBRARY$precMz)
  stopifnot((subset(compSL, precCharge == 2 & isHeavy & grepl("R$",peptide)  )$precMz - subset(SPECTRALLIBRARY, precCharge == 2 & !isHeavy & grepl("R$",peptide))$precMz  )  %>% round  == 5  )
  stopifnot((subset(compSL, precCharge == 4 & isHeavy & grepl("R$",peptide)  )$precMz - subset(SPECTRALLIBRARY, precCharge == 4 & !isHeavy & grepl("R$",peptide))$precMz  )  %>% round  == 3  )
  stopifnot((subset(compSL, precCharge == 3 & !isHeavy & grepl("K$",peptide)  )$precMz - subset(SPECTRALLIBRARY, precCharge == 3 & isHeavy & grepl("K$",peptide))$precMz  )  %>% round  == -3  )

  # all y ion mz should be different
  stopifnot(subset(compSL, ionType =="y")$mz != subset(SPECTRALLIBRARY, ionType =="y")$mz)

  # check fragment mass differences depending on charge state and ion type
  stopifnot((subset(compSL, ionType =="y" & charge == 1 & isHeavy & grepl("R$",peptide)  )$mz - subset(SPECTRALLIBRARY, ionType =="y" & charge == 1 & !isHeavy & grepl("R$",peptide))$mz  )  %>% round  == 10  )
  stopifnot((subset(compSL, ionType =="y" & charge == 2 & isHeavy & grepl("R$",peptide)  )$mz - subset(SPECTRALLIBRARY, ionType =="y" & charge == 2 & !isHeavy & grepl("R$",peptide))$mz  )  %>% round  == 5  )
  stopifnot((subset(compSL, ionType =="y" & charge == 1 & isHeavy & grepl("K$",peptide)  )$mz - subset(SPECTRALLIBRARY, ionType =="y" & charge == 1 & !isHeavy & grepl("K$",peptide))$mz  )  %>% round  == 8  )
  stopifnot((subset(compSL, ionType =="y" & charge == 2 & isHeavy & grepl("K$",peptide)  )$mz - subset(SPECTRALLIBRARY, ionType =="y" & charge == 2 & !isHeavy & grepl("K$",peptide))$mz  )  %>% round  == 4  )

  # all b ion mz should be the same
  stopifnot(subset(compSL, ionType %in% "b")$mz == subset(SPECTRALLIBRARY, ionType =="b")$mz)
  stopifnot(subset(compSL, ionType %in% "b")$mz == subset(SPECTRALLIBRARY, ionType =="b")$mz)

  # # check precMz
  #compSL$peptide = as.character(compSL$peptide)
  #boxplot(compSL$precMz - SPECTRALLIBRARY$precMz ~ paste0(substr(compSL$peptide, nchar(compSL$peptide),nchar(compSL$peptide)),compSL$precCharge, ifelse(compSL$isHeavy,"H","L") ), las=2 )
  #
  # # frag mz
  #boxplot(compSL$mz - SPECTRALLIBRARY$mz ~ paste0(substr(compSL$peptide, nchar(compSL$peptide),nchar(compSL$peptide)),compSL$charge, ifelse(compSL$isHeavy,"H","L"),compSL$ionType  ), las=2 )

  cat("--- testCreateComplementaryIsotopeLibrary:  PASS ALL TEST --- \n")
}


testSpectroDiveExport = function(){

  cat("--- testSpectroDiveExport:  --- \n")
  spectroDiveExport(SPECTRALLIBRARY,TMPXLS)
  cat("--- testSpectroDiveExport:  PASS ALL TEST --- \n")

}

testSpectronautExport = function(){

  cat("--- testSpectronautExport:  --- \n")
  spectronautExport(SPECTRALLIBRARY,TMPXLS)
  cat("--- testSpectronautExport:  PASS ALL TEST --- \n")

}

testSkylineExport = function(){

  cat("--- testSkylineExport:  --- \n")
  spectronautExport(SPECTRALLIBRARY,TMPXLS)

  cat("--- testSkylineExport:  PASS ALL TEST --- \n")

}

testParseTargetsFile = function(){

  cat("--- testParseTargetsFile: --- \n")

  ret = parseTargetsFile(TARGETSFILE)
  stopifnot(length(ret$proteins) == 5 )
  stopifnot(is.na(ret$peptides))

  cat("--- testParseTargetsFile:  PASS ALL TEST --- \n")

}

#run tests

if(T){

  testParseMaxQuantMSMS()
  testCreateLibrarySpectrum()
  #testGetComplementaryIsotopeAnnotatedSpectrum()
  testGetSearchedModifications()
  tesGetIRTModel()
  testGetEmpiricalIRT()
  testGetFragmentSequence()
  #testGetLabelMzShift()
  testGetPeptides()
  testDigestProteome()
  testCreateSpectralLibrary()
  testCreateComplementaryIsotopeLibrary()

  # exports
  testProteotypicPeptideExport()
  testSpectroDiveExport()
  testSpectronautExport()
  testSkylineExport()

  testParseTargetsFile()
}

# GRAPHICS

pdf(TMPPDF)

plotIRTCalibration(IRTMODEL)
barplotPetideCountPerProtein(SPECTRALLIBRARY)
barplotPeptidesPerProtein( cbind(createLibrarySpectrum(TB[7,]),rankingMetric = 1 )[1,] )

cat("CREATED FILE: ", TMPPDF,"\n")


dev.off()







