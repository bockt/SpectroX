
setwd("/Users/ahrnee-adm/dev/R/workspace/SpectroX")
source("R/SpectroX.R")

VERSION <- 2.0

# Use Case I (Find proteotypic peptides):
# Provided a list of proteins
# Select X best peptides per protein from full proteome analysis.
# Peptides are ranked by sum of adjusted fragment intensities
# If not X peptides are found in maxQuant search results, add theoretical peptides (ranked by length) (optional)
# If no protein list is specified use all proteins with at least one peptide
#
#

# 1) parse msms.txt
# 2) apply target protein and/or peptide filters
# 3) calculate iRT times
# 4) create peptide assays
# 5) add theoretical peptides (optional)
# 6) apply peptide sequnece filters
# 7) export

##### PARAMS #########

# input files
MQRESFILE = "./inst/msms.txt"
#MQRESFILE = "~/dev/R/workspace/SpectroUtils/inst/testData/msms.tsv"
#MQRESFILE = "~/dev/R/workspace/SpectroUtils/inst/testData/old/msmsTest.tsv"
#FASTAFILE = "./inst/vibrio_test.fasta"
FASTAFILE = "./inst/uniprot_test.fasta"

# panel filters
TARGETPEPTIDES=NA
TARGETPROTEINS=NA
#TARGETPROTEINS = c("P62258","Q9JII5")

PEPTIDEDPERPROTEN = 5

# peptide filters
PEPTHRS = 0.05
MAXMISCLEAVAGES = 0
PTMREGEXP = "Oxid"
#PTMREGEXP = NA
PEPTIDELENGTHRANGE = c(6,21)
CHARGESTATE = 2:4
#REQUIREDPEPSEQ = "[KR]$"
REQUIREDPEPSEQ = "."
#INVALIDPEPSEQ = NA
#INVALIDPEPSEQ = "(M)|(^[PQ])|([KR]P)"
INVALIDPEPSEQ = "(^[PQ])|([KR]P)"
#LABEL = "light"
#LABEL = "heavy"
LABEL = NA

# assay filters
TRANSITIONS = 5
MINBASEPEAKFRACTION = 0.1
INCLUDENL = F
MINFRAGNB = 3
IONTYPEFILTER = c("b","y")

# misc
KEEPBESTSPECTRUMONLY = T
PROTEASEREGEXP = "[KR](?!P)"

# output files
XLSFILE =  paste0(tempdir(),"/tmp.xls")
PDFFILE = paste0(tempdir(),"/tmp.pdf")

##### PARAMS END #####

### MAIN

# parse maxQuant results
cat("Parsing MaxQuant Results file",MQRESFILE,"\n")
tb = parseMaxQuantMSMS(file = MQRESFILE
                  ,pepThrs = PEPTHRS
                  ,targetPeptides = TARGETPEPTIDES
                  ,targetProteins = TARGETPROTEINS
                  ,filterContaminants = T
                  ,selectedPTMRegExp = PTMREGEXP
                  ,maxMissedCleavages = MAXMISCLEAVAGES
                  ,keepBestSpectrumOnly = T
                  ,chargeState = CHARGESTATE
                  ,label=LABEL
                  )



# peptide sequence filters
# sequence length
keep = tb$Length >= PEPTIDELENGTHRANGE[1] & tb$Length <= PEPTIDELENGTHRANGE[2]
# required peptide sequnece (e.g [KR]$)
keep = keep & grepl(REQUIREDPEPSEQ,tb$Sequence)
# non allowed sequnece fetures e.g (M)|(^[PQ])|([KR]P)
if(!is.na(INVALIDPEPSEQ)) keep =  keep & !grepl(INVALIDPEPSEQ, tb$Sequence)
tb = subset(tb, keep)

# predict iRT
irtModel = getIRTModel(tb)
tb$iRT = getEmpiricalIRT(tb,irtModel$fit)

# create annotated spectra
cat("Parsing Spectra\n")
pbSum <- txtProgressBar(min = 0, max = nrow(tb), style = 3)
spectralLibrary = data.frame()
for(rownb in 1:nrow(tb) ){

  setTxtProgressBar(pbSum, rownb)
  annotSpec = createAnnotatedSpectrum(tb[rownb,])

  ### apply spectrum filters
  # fragment number filter
  filt = with(annotSpec,  (fragmentNb > MINFRAGNB) & grepl(paste(IONTYPEFILTER,collapse="|"), annotSpec$ionType) )
  # neutral loss peak filter
  if(!INCLUDENL) filt = filt &  !annotSpec$isNL
  annotSpec =   subset(annotSpec,filt )

  #fragment intensity as fraction of base peak
  annotSpec$basePeakIntFrac = annotSpec$intensity / suppressWarnings(max(annotSpec$intensity,na.rm=T))
  # filter minimum base peak intensity fraction
  annotSpec = subset(annotSpec,basePeakIntFrac > MINBASEPEAKFRACTION )

  # make sure enough fragments
  if(nrow(annotSpec) > MINFRAGNB){
    # sum of adjusted intensity (basis for ranking)
    annotSpec$adjustedIntensitySum =  annotSpec$adjustedIntensity %>% sum(.,na.rm=T)

    #add to library
    spectralLibrary = rbind(spectralLibrary,annotSpec )
  }
}
# close progress bar
setTxtProgressBar(pbSum, rownb)
close(pbSum)

# parse protein sequence database
if(!is.na(FASTAFILE)){
  proteinDB = read.fasta(FASTAFILE,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
  theoPeptides = digestProteome(proteinDB,peptideLengthRange = PEPTIDELENGTHRANGE
                 ,nbMiscleavages = MAXMISCLEAVAGES
                 ,exclusivePeptides = T
                 ,proteaseRegExp = PROTEASEREGEXP
                 ,trimAC = T)

  # discard peptides already identified
  theoPeptides = subset(theoPeptides, !(peptide %in% spectralLibrary$peptide)  )

  # peptide sequence filters
  # sequence length
  keep = theoPeptides$length >= PEPTIDELENGTHRANGE[1] & theoPeptides$length <= PEPTIDELENGTHRANGE[2]
  # required peptide sequnece (e.g [KR]$)
  keep = keep & grepl(REQUIREDPEPSEQ,theoPeptides$peptide)
  # non allowed sequnece fetures e.g (M)|(^[PQ])|([KR]P)
  if(!is.na(INVALIDPEPSEQ)) keep =  keep & !grepl(INVALIDPEPSEQ,theoPeptides$peptide)
  theoPeptides = subset(theoPeptides, keep)
}

### MAIN END


### EXPORT

# xls
# proteotypic peptide selection
for(protein in TARGETPROTEINS){




}

# graphics
pdf(PDFFILE)
plotIRTCalibration(irtModel)
cat("CREATED FILE: ", PDFFILE,"\n")
dev.off()


### EXPORT END

print(dim(tb))

tb$Proteins

subset(tb, Proteins == "Q9JII5" )

spectralLibrary

print(tb$Modifications)
