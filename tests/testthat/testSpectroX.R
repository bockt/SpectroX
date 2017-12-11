library(tidyverse)
library(magrittr)
library(MASS)
library(seqinr)


setwd("/Users/ahrnee-adm/dev/R/workspace/SpectroX")

source("./R/SpectroX.R")

TESTFILE = "./inst/msms_test.txt"
TB = parseMaxQuantMSMS(TESTFILE)
IRTMODEL = getIRTModel(subset(TB, isIRT))
UPFASTAFILE = "./inst/uniprot_test.fasta"
### read protein db
UPPROTEINDB <- read.fasta(UPFASTAFILE,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
PEPTIDEDF = digestProteome(UPPROTEINDB, peptideLengthRange = c(6,21), dispProgressBar = F)

# unit tests
testParseMaxQuantMSMS = function(){

  cat("--- testParseMaxQuantMSMS: --- \n")

  stopifnot( !(parseMaxQuantMSMS(TESTFILE, label = "light")$Modifications %>% grepl("Arg10|Lys8",.)) )
  stopifnot( sum(parseMaxQuantMSMS(TESTFILE, label = "heavy", keepBestSpectrumOnly=F, selectedPTMRegExp = "Oxid")$Modifications %>% grepl("Arg10|Lys8",.)) == 39 )
  tb = parseMaxQuantMSMS(TESTFILE, pepThrs = 0.01, filterNonExclusivePeptides = T,chargeState = 3, selectedPTMRegExp = "Oxid")
  stopifnot(tb$PEP < 0.01)
  stopifnot(!grepl(";",tb$Proteins))
  #stopifnot(nchar(tb$Sequence)> 15 |  grepl("Biognosys",tb$Proteins) )
  stopifnot(tb$Charge ==  3 |  grepl("Biognosys",tb$Proteins))
  stopifnot( sum( grepl("Oxid",tb$Modifications) ) > 0  )
  stopifnot( sum(rownames(IRTPEPTIDES) %in% tb$Sequence) == 9  )
  # no oxidation
  tb2 = parseMaxQuantMSMS(TESTFILE, pepThrs = 0.01, filterNonExclusivePeptides = F)
  stopifnot( sum( grepl("Oxid",tb2$Modifications) ) == 0  )
  # do not filter non exclusive
  stopifnot(sum(grepl(";",tb2$Proteins)) > 0)
  # no contaminants
  stopifnot( !grepl("CON_",tb2$Proteins) | grepl("Biognosys",tb2$Proteins) )

  tbTargetPep =  parseMaxQuantMSMS(TESTFILE, pepThrs = 0.01, filterNonExclusivePeptides = T, targetPeptides = c("AAVPSGASTGIYEALELR","AFMNNK","AATFGLILDDVSLTHLTFGK") )
  stopifnot(nrow(tbTargetPep) == 1)

  # target peptides
  tbTargetPep =  parseMaxQuantMSMS(TESTFILE, pepThrs = 0.01, filterNonExclusivePeptides = T, targetPeptides = c("AAVPSGASTGIYEALELR","AFMNNK","AATFGLILDDVSLTHLTFGK") )
  stopifnot(nrow(tbTargetPep) == 1)
  # target proteins
  tbTargetProteins =  parseMaxQuantMSMS(TESTFILE, pepThrs = 0.01, filterNonExclusivePeptides = F, targetProteins = c("P67778","P62259","Q9D0K2"), maxMissedCleavages = 2, keepBestSpectrumOnly=F )
  stopifnot(nrow(tbTargetProteins) == 24)

  tbMC = parseMaxQuantMSMS(TESTFILE, pepThrs = 1, maxMissedCleavages = 2 )

  nrow(tbMC)

  # keepBestSpectrumOnly
  stopifnot(table(paste(tbMC$Sequence,tbMC$Modifications)) %>% max ==1)
  stopifnot((parseMaxQuantMSMS(TESTFILE, pepThrs = 1, maxMissedCleavages = 2, keepBestSpectrumOnly = F ) %>% nrow) > nrow(tbMC))

  stopifnot( tbMC$`Missed cleavages` %>% max == 2)
  cat("--- testParseMaxQuantMSMS: PASS ALL TEST --- \n")

}

testGetSearchedModifications = function(){
  cat("--- testGetSearchedModifications: --- \n")

  #stopifnot( c("Arg10","Lys8","Oxidation (M)") %in%  getSearchedModifications(TB))
  stopifnot( c("Oxidation (M)") %in%  getSearchedModifications(TB))
  cat("--- testGetSearchedModifications: PASS ALL TEST --- \n")
}

testCreateAnnotatedSpectrum <- function(){

  cat("--- testParseFragmentInformation: --- \n")
  annotSpec = createAnnotatedSpectrum(TB[7,])
  stopifnot(dim(annotSpec) == c(17,18) )
  stopifnot( c("a","b","y") == levels(annotSpec$ionType) )
  stopifnot(annotSpec$adjustedIntensity > annotSpec$intensity    )
  stopifnot(annotSpec$precCharge == TB[7,]$Charge)
  cat("--- testParseFragmentInformation: PASS ALL TEST --- \n")
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
  stopifnot(!is.na(createAnnotatedSpectrum(tb[9,])$iRT) )
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

testGetLabelMzShift <- function(){

  cat("--- testGetLabelMzShift: --- \n")

  ### ONLY LABELLED AT THE TERMINI

  stopifnot(getLabelMzShift(aaSeq = "PETIDEK", charge=1, isHeavy=T) == -8.014199)
  stopifnot(getLabelMzShift(aaSeq = "PEPTIDER", charge=1, isHeavy=T) == -10.008269)
  stopifnot(getLabelMzShift(aaSeq = "PEKPTIDEK", charge=1, isHeavy=F) == 8.014199)
  stopifnot(getLabelMzShift(aaSeq = "PERPTIDER", charge=2, isHeavy=F) == 10.008269/2)

  cat("--- testGetLabelMzShift:  PASS ALL TEST --- \n")

}

testGetComplementaryIsotopeAnnotatedSpectrum <- function(){

  cat("--- testGetComplementaryIsotopeAnnotatedSpectrum: --- \n")

  annotSpec = createAnnotatedSpectrum(TB[9,])
  compAnnotSpec = getComplementaryLabelSpectrum(annotSpec)

  stopifnot(subset(annotSpec, ionType == "y")$mz < subset(compAnnotSpec, ionType == "y")$mz)
  stopifnot(subset(annotSpec, ionType == "b")$mz == subset(compAnnotSpec, ionType == "b")$mz)
  stopifnot(subset(annotSpec, ionType == "a")$mz == subset(compAnnotSpec, ionType == "a")$mz)
  stopifnot(compAnnotSpec$isHeavy == !annotSpec$isHeavy)

  cat("--- testGetComplementaryIsotopeAnnotatedSpectrum:  PASS ALL TEST --- \n")

}

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

#run tests

if(T){

  testParseMaxQuantMSMS()
  testCreateAnnotatedSpectrum()
  testGetComplementaryIsotopeAnnotatedSpectrum()
  testGetSearchedModifications()
  tesGetIRTModel()
  testGetEmpiricalIRT()
  testGetFragmentSequence()
  testGetLabelMzShift()
  testGetPeptides()
  testDigestProteome()
}


# GRAPHICS
tmpPDF = paste0(tempdir(),"/tmp.pdf")
pdf(tmpPDF)
plotIRTCalibration(IRTMODEL)
cat("CREATED FILE: ", tmpPDF)
dev.off()





























