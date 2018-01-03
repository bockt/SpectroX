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

# predict iRT
irtModel = getIRTModel(tb)
tb$iRT = getEmpiricalIRT(tb,irtModel$fit)

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
pdf(uO$PDFFILE)
parDefault = par()
plotIRTCalibration(irtModel)

# barplot peptide count per protien
barplotPetideCountPerProtein(spectralLibrary)
#  plot adj. intensity vs peptide (per protein)
par(mfrow=c(2,2), mar=c(5,5,5,5))

X = spectralLibrary[c("protein","peptide","ptm","rankingMetric","precCharge")] %>%
  split(.$protein) %>%
  map(~  barplotPeptidesPerProtein(.x %>% unique , ptmRegExp = uO$PTMREGEXP, rankingMetric=uO$RANKINGMETRIC))


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

