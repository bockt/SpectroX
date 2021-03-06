#!/usr/bin/Rscript

#suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(seqinr))

# set this
PATHTOSPECTROX = "/Users/ahrnee-adm/dev/R/workspace/SpectroX"
if(!file.exists(PATHTOSPECTROX)){
  PATHTOSPECTROX = "/home/pcfuser/R/SpectroX"
}

#setwd(PATHTOSPECTROX)
#source("R/SpectroX.R")
#source("R/UserOptions.R")

source(paste0(PATHTOSPECTROX,"/R/SpectroX.R"))
source(paste0(PATHTOSPECTROX,"/R/UserOptions.R"))

VERSION <- 2.0

# Use Case I (Find proteotypic peptides):
# Provided a list of proteins
# Select X best peptides per protein from full proteome analysis.
# Rank peptides by sum of adjusted fragment intensities
# If not X peptides are found in maxQuant search results, add theoretical peptides (ranked by length) (optional)
# If no protein list is specified use all proteins with at least one peptide
#

# Use Case II:
# Provided a list a peptides or proteins
# Rank peptides by sum of adjusted fragment intensities
# Optional: created complementary isotope spectra
# Export Spectral Library (panel) in spectrodive, spectronaut or skyline format

# 1) parse msms.txt
# 2) apply target protein and/or peptide filters
# 3) calculate iRT times
# 4) create peptide assays
# 5) add theoretical peptides (optional)
# 6) apply peptide sequnece filters
# 7) export

##### PARAMS #########


### USER CMD LINE OPTIONS
cmdlineOptions <- getCMDLineOptions(version=VERSION)
uO = getUserOptions(cmdlineOptions)
### USER CMD LINE OPTIONS END

# panel filters

#### PARSE PEPPROTLISTFILE

##### PARAMS END #####

### MAIN

# parse maxQuant results
cat("Parsing MaxQuant Results file",uO$MQRESFILE,"\n")
# apply table filters
tb = parseMaxQuantMSMS(file = uO$MQRESFILE
                  ,pepCutoff = uO$PEPCUTOFF
                  ,targetPeptides = uO$TARGETPEPTIDES
                  ,targetProteins = uO$TARGETPROTEINS
                  ,filterContaminants = T
                  ,selectedPTMRegExp = uO$PTMREGEXP
                  ,maxMissedCleavages = uO$MAXMISCLEAVAGES
                  ,keepBestSpectrumOnly = T
                  ,chargeState = uO$CHARGESTATE
                  ,label=uO$LABEL
                  ,pepLength =uO$PEPTIDELENGTHRANGE
                  )

keep = rep(T,nrow(tb))
# apply peptide sequence filters unless target peptide list provided
if(is.na(uO$TARGETPEPTIDES)[1]){
  keep = grepl(uO$REQUIREDPEPSEQ,tb$Sequence)
  # non allowed sequnece fetures e.g (M)|(^[PQ])|([KR]P)
  if(!is.na(uO$INVALIDPEPSEQ)) keep =  keep & !grepl(uO$INVALIDPEPSEQ, tb$Sequence)
  tb = subset(tb, keep | isIRT)
}
# stop if all psms are filtered out
if(nrow(tb) == 0 )stop("0 PSMs passing filtering criteria")

# predict iRT
irtModel = getIRTModel(tb)
if(is.null(irtModel$fit)){ # not enough iRT
  tb$iRT = rep(NA,nrow(tb))
}else{
  tb$iRT = getEmpiricalIRT(tb,irtModel$fit)
}

if(uO$VERBOSE){
  cat("NB PSMs: ",nrow(tb),"\n")
}


# create spectral library
spectralLibrary = createSpectralLibrary(tb
                      , minFragNb= uO$MINFRAGNB
                      , minNbTransitions= min(uO$NBTRANSITIONSRANGE)
                      , maxNbTransitions =  max(uO$NBTRANSITIONSRANGE)
                      , minBasePeakFraction = uO$MINBASEPEAKFRACTION
                      , ionTypeFilter = uO$IONTYPEFILTER
                      , includeNL = uO$INCLUDENL
                      ,rankingMetric = uO$RANKINGMETRIC
                      )

# parse protein sequence database
if(!is.na(uO$FASTAFILE)){
  proteinDB = read.fasta(uO$FASTAFILE,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
  theoPeptides = digestProteome(proteinDB,peptideLengthRange = uO$PEPTIDELENGTHRANGE
                 ,nbMiscleavages = uO$MAXMISCLEAVAGES
                 ,exclusivePeptides = T
                 ,proteaseRegExp = uO$PROTEASEREGEXP
                 ,trimAC = T)

   # discard peptides already identified
  theoPeptides = subset(theoPeptides, !(peptide %in% spectralLibrary$peptide)  )

  # peptide sequence filters
  # sequence length
  keep = theoPeptides$length >= min(uO$PEPTIDELENGTHRANGE) & theoPeptides$length <= max(uO$PEPTIDELENGTHRANGE)
  # required peptide sequnece (e.g [KR]$)
  keep = keep & grepl(uO$REQUIREDPEPSEQ,theoPeptides$peptide)
  # non allowed sequnece fetures e.g (M)|(^[PQ])|([KR]P)
  if(!is.na(uO$INVALIDPEPSEQ)) keep =  keep & !grepl(uO$INVALIDPEPSEQ,theoPeptides$peptide)
  theoPeptides = subset(theoPeptides, keep)
}

# select top peptide variants (peptide+ ptm ) per protein
selPeptideVariants = spectralLibrary[c("protein","peptide","ptm","rankingMetric","isIRT")]%>%
   unique %>% group_by(protein,peptide,ptm, isIRT) %>% top_n(1,rankingMetric) %>%
   split(.$protein) %>%
   map_df(~ data.frame(peptide=.x$peptide
                       , ptm = .x$ptm
                       , protein = .x$protein
                       , isIRT = .x$isIRT
                       , rankingMetric= .x$rankingMetric
                       , rank = length(.x$rankingMetric) - rank(.x$rankingMetric)+1
   ))



selPeptideVariants = subset(selPeptideVariants,rank <= uO$PEPTIDEDPERPROTEN | isIRT )
# apply peptide ptm filter
spectralLibrary = spectralLibrary[ paste0(spectralLibrary$peptide,spectralLibrary$ptm) %in% paste0(selPeptideVariants$peptide,selPeptideVariants$ptm),]

### MAIN END

### EXPORT

# graphics
pdf(uO$PDFFILE )
parDefault = par()
par(cex.axis = 1.3, cex.lab =1.3, mfrow=c(2,2))
if(!is.null(irtModel$fit)) plotIRTCalibration(irtModel)

# barplot peptide count per protein
barplotPetideCountPerProtein(spectralLibrary)

# add ionTypeTag e.g. yNL++
spectralLibrary$ionTypeTag =  with(spectralLibrary,  paste0(ionType
                                  ,ifelse(isNL,"NL","")
                                  ,lapply(spectralLibrary$charge, function(t){paste(rep("+",t), collapse="")} ) %>% unlist
))

# # fragmens per spectrum
barplot((group_by(spectralLibrary,mqResIdx) %>% dplyr::count())$n %>% table
         , col = "blue"
         , ylab ="PSM Count"
         #, cex.axis = 1.5
         , cex.names =1.25
         #, cex.lab =1.5
         , xlab  = "Fragments per Assay"
       )

# charge state distrib
psmPerCharge = group_by(spectralLibrary, precCharge)  %>% summarise(nbPSM = length(unique(mqResIdx)) )
barplot(psmPerCharge$nbPSM
        ,names = psmPerCharge$precCharge
        , col = "blue"
        , ylab ="PSM Count"
        #, cex.axis = 1.5
        , cex.names =1.25
        #, cex.lab =1.5
        ,xlab  = "Precursor Charge"
)

# plot fragment type distribution
# precursor charge states with at least 5 psms
#charges = (group_by(spectralLibrary,precCharge) %>% summarise(count = length(precCharge) ) %>% filter(count > 10 ))$precCharge
#charges = subset(psmPerCharge, nbPSM > 5 )$precCharge
charges = psmPerCharge$precCharge

par(mfrow = c( max(1,ceiling(length(charges) / 2)),2) )
X = subset(spectralLibrary, precCharge %in% charges ) %>%
  split(.$precCharge) %>%
  # map(~  barplotPeptidesPerProtein(.x %>% unique , ptmRegExp = uO$PTMREGEXP, rankingMetric=uO$RANKINGMETRIC))
  map(~  barplot(100* table(.x$ionTypeTag) / nrow(.x), main = paste0("Precursor Charge ",.x$precCharge[1],"+" )
                 , col = "blue"
                 , ylab ="%"
                 #, cex.axis = 1.5
                 , cex.names =1.25
                 #, cex.lab =1.5
                 ,las=2
          )
  )
par(mfrow = c(1,1), cex.lab = 1, cex.axis =1 )

#  plot adj. intensity vs peptide (per protein)
par(mfrow=c(2,2), mar=c(5,5,5,5))

# plot ptm color legend
ptmCol = getPTMColors(levels(spectralLibrary$ptm) )
plot(0,0,type="n", xaxt="n", yaxt="n", axes=FALSE,xlab="",ylab="")
legend("center",fill =ptmCol$col, legend=rownames(ptmCol),cex=0.7, bty="n")


X = spectralLibrary[c("protein","peptide","ptm","rankingMetric","precCharge")] %>%
  split(.$protein) %>%
  map(~  barplotPeptidesPerProtein(.x %>% unique
                                   , rankingMetric=uO$RANKINGMETRIC
                                   ,ptmCol=ptmCol))


cat("CREATED FILE: ", uO$PDFFILE,"\n")
graphics.off()

# xls
if(uO$PPEXPORT ){

  if(!exists("theoPeptides")){
    theoPeptides = NA
    #warning("No predicted peptides added")
  }

  # proteotypic peptide selection
  proteotypicPeptideExport(spectralLibrary = spectralLibrary
                           ,targetProteins = uO$TARGETPROTEINS
                           , nbPeptidesPerProtein = uO$PEPTIDEDPERPROTEN
                           , theoPeptides = theoPeptides
                           , outFile =uO$XLSFILE
  )
}

if(uO$COMPLEMENTARYISOTOPEASSAYS){
  spectralLibrary = rbind(spectralLibrary, subset(createComplementaryIsotopeLibrary(spectralLibrary), !(isIRT & isHeavy)  ))
  # group light heavy
  spectralLibrary =arrange(rbind(spectralLibrary,createComplementaryIsotopeLibrary(spectralLibrary)),mqResIdx,precMz)
}

if(uO$SPECTRODIVEEXPORT){
  spectroDiveExport(spectralLibrary,file = uO$XLSFILE )
}

if(uO$SPECTRONAUTEXPORT){
  spectronautExport(spectralLibrary,file = uO$XLSFILE )
}

if(uO$SKYLINEEXPORT){
  skylineExport(spectralLibrary,file = uO$XLSFILE )
}

### EXPORT END

