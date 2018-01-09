# TODO: Add comment
#
# Author: erikahrne
###############################################################################

suppressPackageStartupMessages(library(optparse))



#' Parse range in format 2:4 and return integer vector. Inthis case 2,3,4
#' @param rangeStr e.g. '2:4'
#' @return vectot of integers
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
parseRange =  function(rangeStr){

  if(!grepl("([1-9]{1,2})(\\:){0,1}([1-9]{0,2})",rangeStr) ) stop("Invalid range: ", rangeStr)
  r = eval(parse(text = rangeStr))
  if(!(is.integer(r) | is.numeric(r)) ) stop("Invalid range: ", rangeStr)
  return(c(min(r),max(r)))
}

#' Define and User Specified Command Line Options
#' @param version SpectroX version number
#' @return list() cmd line options
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getCMDLineOptions = function(version=version){

  ### CMD OPTIONS
  option_list <- list(

    ### I/O
    make_option(c("-i", "--inputFile"), type="character", default=NA,
                help="I/O: MaxQuant msms.txt input file (REQUIRED)"
    ),
    make_option(c("-p", "--pepProtListFile"), type="character", default=NA,
                help="I/O: comma-separated file listing target peptides/proteins acs (use header 'Peptide' or 'Protein')
                If no file is provided, all identified peptides will be used."
    ),
    make_option(c("-f", "--fastaFile"), type="character", default=NA,
                help="I/O: protein database .fasta file (only relevant for 'Proteotypic peptide selection' workflow)"
    ),
    make_option(c("-o", "--outputDir"), type="character", default="./",
                help="I/O:  Results Output Directory [default %default]"
    ),

    make_option(c("-l", "--resultsFileLabel"), type="character", default="spectroXInput",
                help="I/O: .pdf & .xls file labels (prefix) [default %default]"
    ),

    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="I/O: Print extra output [default %default]"
    ),

    ### I/O END

    # PEPTIDE (-P)

    make_option(c("--PPeptidesPerProtein"), type="integer", default=5,
                help="PEPTIDE: --PP Targeted number of peptides per protein [default %default]"
    ),

    make_option(c("--PModificationRegEx"), type="character", default=NA,
                help="PEPTIDE: --PM Allowed PTMs. E.g. \"Oxid|Ac\" -> Oxidation or Acetylation [default %default]"
    ),

    make_option(c("--PEpCutoff"), type="double", default=0.05,
                help="PEPTIDE: --PE Posterior Error Probability Cutoff.  [0-1] [default %default]"
    ),

    make_option(c("--PNumberMisCleavages"), type="integer", default=0,
                help="PEPTIDE: --PN Maximum number of missed cleavages. [default %default]"
    ),

    make_option(c("--PCleavageEnzymeRegEx"), type="character", default="[KR](?!P)",
                help="PEPTIDE: --PC Protease cleavage specificity regular expression. [default %default]"
    ),

    make_option(c("--PLengthRange"), type="character", default="5:30",
                help="PEPTIDE: --PL Peptide length filter. Example '5:30' keeps peptides with [5,30] amino acids  . [default no filter]"
    ),

    make_option(c("--PZRange"), type="character", default=NA,
                help="PEPTIDE: --PZ Peptide charge state filter. Example '2:4' keeps charge states 2 ,3 and 4  . [default no filter]"
    ),

    make_option(c("--PRequiredSeqRegEx"), type="character", default="[KR]$",
                help="PEPTIDE: --PR Required Peptide Sequence. E.g. '[KR]$' -> tryptic c-term [default %default]"
    ),

    make_option(c("--PInvalidSeqRegEx"), type="character", default="(M)|(^[EQ])|([KR]P)",
                help="PEPTIDE: --PC Invalid Peptide Sequence.
                Example '(M)|(^[EQ])|([KR]P)' ->
                No Met containing peptides.
                No Peptides starting to with Glu or Gln, prone to pyro-Glu.
                No Peptides having subsequence Lys or Arg followed by Pro
                [default %default]"
    ),

    make_option(c("--PHeavyIsotopeLabel"), action="store_true", default=NA,
                help="PEPTIDE: --PH Keep Heavy Lys and Arg peptides only.
                [default no filter]
                --PH -> heavy only
                --PH F -> light only"
    ),

    # PEPTIDE (-P) END

    # SPECTRUM (--S)

    make_option(c("--STransitionNbRange"),type="character", default=5,
                help="SPECTRUM: --SF Spectrum transition number range.
                Example '5:6' .
                Min 5 transitions
                Max 6 transitions
                [default %default]"
    ),

    make_option(c("--SBasePeakFraction"), type="double", default=0.0,
                help="SPECTRUM: --SB Filter out peaks with intensity < set value.  [0-1] [default %default]"
    ),

    make_option(c("--SNeutralLossFragments"), action="store_true", default=FALSE,
                help="SPECTRUM: --SN Include Neutral Loss Fragments  [default %default]"
    ),

    make_option(c("--SFragmentAnnotationNumberMin"), type="integer", default=3,
                help="SPECTRUM: --SF Minimum fragment number
                3 -> e.g exlude  a1,a2,b1,b2,x1,x2,y1,y2
                [default %default]"
    ),

    make_option(c("--SIonTypesAllowed"),type="character", default="abxy",
                help="SPECTRUM: --SI allows ion types.
                Example: 'yb' -> B and Y ions only
                [default %default]"
    ),

    make_option(c("--SComplementaryIsotopeSpectrum"), action="store_true", default=TRUE,
                help="SPECTRUM: --SC create complementary isotope spectra.
                I.e. if Heavy Lys peptide was identified, also export light peptide assay.
                --SC F -> flag deactivated i.e. light or heavy only
                [default %default]"
    ),

    # SPECTRUM (--S) END

    # RANKING (--R)

    make_option(c("--RRankingMetric"), type="integer", default=1,
                help="RANKING: --RR [default %default]
                1) adjustedIntensitySum
                2) precApexIntensity
                3) psmScore"
    ),
    # RANKING (--R) END


    #EXPORT (--E)
    make_option(c("--EXportFormats"), type="integer", default=1,
                help="EXPORT: --EX [default %default]
                1) SpectroDive compatible
                2) Spectronaut
                3) Skyline
                4) Proteotypic peptide selection"

    )
    #EXPROT (--E) END

    )

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser(prog=paste("spectroX",version), option_list=option_list))

	printOptions(cmdOpt, header="CMD LINE Options" )

	### CMD OPTIONS END

# TEST END
	return(cmdOpt)

}

printOptions = function(cmdOpt, header="Options"){
   cat("\n-------- SpectroX ",header," ------------\n")
   for(opt in names(cmdOpt)){
     cat(opt," ",cmdOpt[[opt]],"\n")
   }
   cat("--------------------------------------\n")

}

#' check and return User Specified Command Line Options
#' @param list()  command line options
#' @return list() user options
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getUserOptions = function(cmdlineOptions){

  uO = list()

  ############################## input files ##############################
  uO$MQRESFILE = cmdlineOptions$inputFile
  if(!file.exists(uO$MQRESFILE)) stop("File not found: ",uO$MQRESFILE,"\n")

  uO$FASTAFILE = cmdlineOptions$fastaFile
  if(!is.na(uO$FASTAFILE) & !file.exists(uO$FASTAFILE)) stop("File not found: ",uO$FASTAFILE,"\n")

  uO$OUTDIR = cmdlineOptions$outputDir
  if(is.na(uO$OUTDIR)) uO$OUTDIR  = dirname( uO$MQRESFILE)
  if(!file.exists(uO$OUTDIR)) stop("outputDir not found: ",uO$OUTDIR,"\n")
  #if(!is.na(uO$OUTDIR) & !file.exists(uO$OUTDIR)) stop("outputDir not found: ",uO$OUTDIR,"\n")

  uO$PEPPROTLISTFILE = cmdlineOptions$pepProtListFile
  if(!is.na(uO$PEPPROTLISTFILE) & !file.exists(uO$PEPPROTLISTFILE)) stop("File not found: ",uO$PEPPROTLISTFILE,"\n")

  ##############################  peptide filters ##############################
  uO$PEPTIDEDPERPROTEN = cmdlineOptions$PPeptidesPerProtein
  if(uO$PEPTIDEDPERPROTEN < 1) stop("Invalid PPeptidesPerProtein: ",uO$PEPTIDEDPERPROTEN,"\n")

  uO$PTMREGEXP = cmdlineOptions$PModificationRegEx

  uO$PEPCUTOFF = cmdlineOptions$PEpCutoff
  if(uO$PEPCUTOFF < 0 | uO$PEPCUTOFF > 1 ) stop("Invalid PEpCutoff: ",uO$PEPCUTOFF,"\n")

  uO$MAXMISCLEAVAGES = cmdlineOptions$PNumberMisCleavages
  if(uO$MAXMISCLEAVAGES < 0) stop("Invalid PNumberMisCleavages: ",uO$MAXMISCLEAVAGES,"\n")

  uO$PROTEASEREGEXP = cmdlineOptions$PCleavageEnzymeRegEx

  if(is.na(cmdlineOptions$PLengthRange)){
    uO$PEPTIDELENGTHRANGE = c(0,Inf)
  }else{
    uO$PEPTIDELENGTHRANGE  = parseRange(cmdlineOptions$PLengthRange)
  }

  if(is.na(cmdlineOptions$PZRange)){
    uO$CHARGESTATE = c(0,Inf)
  }else{
    uO$CHARGESTATE  = parseRange(cmdlineOptions$PZRange)
  }

  uO$REQUIREDPEPSEQ = cmdlineOptions$PRequiredSeqRegEx
  uO$INVALIDPEPSEQ = cmdlineOptions$PInvalidSeqRegEx

  uO$LABEL = cmdlineOptions$PHeavyIsotopeLabel
  if(!is.na(uO$LABEL)) uO$LABEL = ifelse(uO$LABEL,"heavy","light")

  ############################## spectrum filters ##############################
  if(is.na(cmdlineOptions$STransitionNbRange)){
    uO$NBTRANSITIONSRANGE = c(0,Inf)
  }else{
    uO$NBTRANSITIONSRANGE  = parseRange(cmdlineOptions$STransitionNbRange)
  }

  uO$MINBASEPEAKFRACTION = cmdlineOptions$SBasePeakFraction
  if(uO$MINBASEPEAKFRACTION < 0 | uO$MINBASEPEAKFRACTION > 1 ) stop("Invalid SBasePeakFraction: ",uO$MINBASEPEAKFRACTION,"\n")

  uO$INCLUDENL = cmdlineOptions$SNeutralLossFragments

  uO$MINFRAGNB = cmdlineOptions$SFragmentAnnotationNumberMin
  if(uO$MINFRAGNB < 1 ) stop("Invalid SFragmentAnnotationNumberMin: ",uO$MINFRAGNB,"\n")


  uO$IONTYPEFILTER = c("a","b","c","x","y","z")
  if(!is.na(cmdlineOptions$SIonTypesAllowed)){
    itSel = (lapply(strsplit(cmdlineOptions$SIonTypesAllowed %>% tolower() ,"")[[1]],function(t){grepl(t,uO$IONTYPEFILTER)}) %>% data.frame %>% rowSums) >0
    uO$IONTYPEFILTER = uO$IONTYPEFILTER[itSel]
  }
  if(length(uO$IONTYPEFILTER) == 0)  stop("Invalid SIonTypesAllowed: ",cmdlineOptions$SIonTypesAllowed,"\n")

  uO$COMPLEMENTARYISOTOPEASSAYS = cmdlineOptions$SComplementaryIsotopeSpectrum

  ############################## peptide or protein  targets ##############################
  uO$TARGETPEPTIDES=  NA
  uO$TARGETPROTEINS= NA
  if(!is.na(uO$PEPPROTLISTFILE)){
    targets = parseTargetsFile(uO$PEPPROTLISTFILE)
    uO$TARGETPEPTIDES=  targets$peptides
    uO$TARGETPROTEINS= targets$proteins
    # set uO$PEPTIDEDPERPROTEN to Inf if target peptide list provided
    if(!is.na(uO$TARGETPEPTIDES[1])) uO$PEPTIDEDPERPROTEN = Inf
  }

  ############################## raning ##############################

  # TEMP
  #uO$RANKINGMETRIC = "adjustedIntensitySum"
  #uO$RANKINGMETRIC = "psmScore"
  #uO$RANKINGMETRIC = "RANKINGMETRIC"

  uO$RANKINGMETRIC = "adjustedIntensitySum"
  if(cmdlineOptions$RRankingMetric == 2) uO$RANKINGMETRIC = "precApexIntensity"
  if(cmdlineOptions$RRankingMetric == 3) uO$RANKINGMETRIC = "psmScore"
  if(cmdlineOptions$RRankingMetric > 3)  stop("Invalid RRankingMetric: ",cmdlineOptions$RRankingMetric,"\n")

  ############################## output files ##############################

  uO$PPEXPORT = F
  uO$SPECTRODIVEEXPORT = F
  uO$SPECTRONAUTEXPORT = F
  uO$SKYLINEEXPORT = F
  if(cmdlineOptions$EXportFormats == 4) uO$PPEXPORT = T
  if(cmdlineOptions$EXportFormats == 1) uO$SPECTRODIVEEXPORT = T
  if(cmdlineOptions$EXportFormats == 2) uO$SPECTRONAUTEXPORT = T
  if(cmdlineOptions$EXportFormats == 3) uO$SKYLINEEXPORT = T
  if(sum(uO[grepl("EXPORT$", names(uO))] %>% unlist ) == 0 )stop("Invalid EXportFormats: ",cmdlineOptions$EXportFormats,"\n")


  uO$XLSFILE =  paste0(uO$OUTDIR,"/", cmdlineOptions$resultsFileLabel,".xls")
  uO$PDFFILE =  paste0(uO$OUTDIR,"/", cmdlineOptions$resultsFileLabel,".pdf")

  uO$VERBOSE = cmdlineOptions$verbose

  if(uO$VERBOSE) printOptions(uO)

  return(uO)

}
