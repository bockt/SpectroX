#
# erik.ahrne@unibas.ch
#

# IRTPEPTIDES = c("LGGNEQVTR"
#                 ,"GAGSSEPVTGLDAK"
#                 ,"VEATFGVDESNAK"
#                 ,"YILAGVENSK"
#                 ,"TPVISGGPYEYR"
#                 ,"TPVITGAPYEYR"
#                 ,"DGLDAASYYAPVR"
#                 ,"ADVTPADFSEWSK"
#                 ,"GTFIIDPGGVIR"
#                 ,"GTFIIDPAAVIR"
#                 ,"LFLQFGAQGSPFLK")

#' @export
IRTPEPTIDES <- data.frame(refIRT = sort(c(0.00,0.00
                                     ,54.62, 54.62
                                     ,100.00, 100.00
                                     ,42.26, 42.26
                                     ,87.23, 87.23
                                     ,-24.92, -24.92
                                     ,70.52, 70.52
                                     ,28.71,  28.71
                                     ,33.38, 33.38
                                     ,12.39, 12.39
                                     ,19.79,19.79))
                     ,row.names=c("LGGNEQVTR", "LGGNETQVR"
                                  ,"GAGSSEPVTGLDAK", "AGGSSEPVTGLADK"
                                  ,"VEATFGVDESNAK", "VEATFGVDESANK"
                                  ,"YILAGVENSK", "YILAGVESNK"
                                  ,"TPVISGGPYEYR","TPVISGGPYYER"
                                  ,"TPVITGAPYEYR", "TPVITGAPYYER"
                                  ,"DGLDAASYYAPVR", "GDLDAASYYAPVR"
                                  ,"ADVTPADFSEWSK", "DAVTPADFSEWSK"
                                  ,"GTFIIDPGGVIR", "TGFIIDPGGVIR"
                                  ,"GTFIIDPAAVIR", "GTFIIDPAAIVR"
                                  ,"LFLQFGAQGSPFLK", "FLLQFGAQGSPLFK"))


#' Parse MaxQuant msms.txt
#' @param file path
#' @param pepThrs numeric default 0.05
#' @param targetPeptides character default NA
#' @param targetProteins character default NA
#' @param filterContaminants TRUE
#' @param contaminantRegExp '^CON_'
#' @param selectedPTMRegExp NA
#' @param filterNonExclusivePeptides default TRUE
#' @param chargeState default [1,10]
#' @param label (Arg10 and Lys8 only filter) options NA (deafault,no filter),  'light','heavy'
#' @param keepBestSpectrumOnly  keep top-scoring spectrum only, per peptidoform, default TRUE
#' @param requiredSequenceRegExp dafualt '.' (no filter)
#' @return data.frame of maxQuant psm level results
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
parseMaxQuantMSMS = function(file
                             , pepThrs=0.05
                             , targetPeptides=NA
                             , targetProteins=NA
                             , filterContaminants = T
                             , contaminantRegExp = '^CON_'
                             , selectedPTMRegExp = NA
                             , filterNonExclusivePeptides = T
                             , chargeState = 1:10
                             , label = NA
                             , maxMissedCleavages = 0
                             , keepBestSpectrumOnly = T
                             ,requiredPepSeqRegExp = "."

){

  tb = suppressMessages(read_tsv(file))

  # is irt peptide
  tb$isIRT = paste(IRTPEPTIDES %>% rownames , collapse="|") %>% grepl(.,tb$`Modified sequence`)

  # target peptide filter
  if(!is.na(targetPeptides[1])){
    # targetPeptides = c("AAAAAAAAALGVR", "AAAAADEDEEADEEDEAEAGGMGPQAR" )
    tb = subset(tb, Sequence %in% targetPeptides  )
  }
  #Filter out non exclusive peptides
  if(filterNonExclusivePeptides){
    tb = subset(tb, !grepl("\\;",Proteins))
  }else if(!is.na(targetPeptides[1])){ # warn if target peptides are not exclusive
    nonExclPeptides =  subset(tb, grepl("\\;",Proteins))
    if(nrow(nonExclPeptides) > 0 ){
      cat("Warning - Found Non-Exclusive Peptides:")
      nonExclPeptides[c("Sequence","Proteins")] %>% as.data.frame() %>% print
    }
  }

  # non peptide filters
  keep = rep(T, nrow(tb))
  # target protein filter
  if(!is.na(targetProteins[1])){
    #targetProteins = c("P46099", "Q99J95" )
    targetRegex = paste(targetProteins,collapse = "|")
    keep = keep & grepl(targetRegex , tb$Proteins)
  }
  # CON filter
  if(filterContaminants) keep = keep & !grepl("CON\\_",tb$Proteins)
  # pep thrs
  keep = keep & (tb$PEP <  pepThrs)
  # decoy filter
  keep = keep & (is.na(tb$Reverse) | tb$Reverse != "+" )
  # min peptide length
  # keep = keep &  nchar(tb$Sequence) > minPepLength
  # keep charge state
  keep = keep &  (tb$Charge %in% chargeState)
  # keep modifs
  #selectedPTMRegExp = "Oxidation"
  # if(!is.na(selectedPTMRegExp)){
  #   keep =  keep & (paste(c("Unmodified","Arg10","Lys8",selectedPTMRegExp),collapse="|") %>% grepl(.,tb$Modifications))
  # }else{
  #   keep =  keep & (paste(c("Unmodified","Arg10","Lys8"),collapse="|") %>% grepl(.,tb$Modifications))
  # }

  # get all modifications seached for
  excludedPTM = getSearchedModifications(tb)
  cat("INFO: allowed PTMs:", ifelse(is.na(excludedPTM[grepl(selectedPTMRegExp,excludedPTM)]),"None", paste(excludedPTM[grepl(selectedPTMRegExp,excludedPTM)],collapse=":")  ),"\n")
  # e.g. "Oxidation (M)". Is "Oxidation (M)" allowed. Check if it matches regexp
  if(!is.na(selectedPTMRegExp)) excludedPTM = excludedPTM[!grepl(selectedPTMRegExp,excludedPTM)]

  # if excludedPTM vector contains any mods. only keep peptides not affected by these ptms
  if(length(excludedPTM) > 0) keep =  keep & rowSums(tb[excludedPTM]) == 0

  # light heavy filter
  if(!is.na(label)){
    isHeavyIsotope = grepl("Arg10|Lys8",tb$Modifications)
    if(label == "light"){ # keep light only
      keep =  keep & !isHeavyIsotope
    }else{ # keep heavy only
      keep =  keep & isHeavyIsotope
    }
  }

  # missed cleavages
  keep =  keep & tb$`Missed cleavages` <= maxMissedCleavages

  # DO NOT discard irts
  #keep = keep | ((paste(IRTPEPTIDES, collapse="|") %>% grepl(.,tb$`Modified sequence`)) &  (tb$PEP <  pepThrs))
  keep = keep | (tb$isIRT & (tb$PEP <  pepThrs))

  # apply filter
  tb = subset(tb, keep)

  # filterBestSpectrum
  # keep top-scoring spectrum only, per peptidoform
  if(keepBestSpectrumOnly) tb = group_by(tb, Sequence, Modifications) %>% filter(rank(Score,ties.method = "random") == length(Score))

  # casting
  suppressWarnings(tb$`Precursor Apex Fraction` %<>% as.numeric)

  return(tb)

}

#' Get list of variable modifiections condsidered by MaxQuant
#' @param df tible or data.frame
#' @return character vector of modification names (does not return "Arg10","Lys8" labels )
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getSearchedModifications = function(tb){

  mods = names(tb)[grepl("Probabilities$",names(tb))] %>% gsub(" Probabilities$","",.)
  # remove "Arg10"         "Lys8"
  mods = mods[!(mods %in% c("Arg10","Lys8"))]
  return(mods)
}

#' Parse spectrum framgment match information listed in MaxQuant 'Identifications' results
#' @param psm a row of MaxQuant 'Identifications' data.frame
#' @return list with  attributes intensity, ionType, charge, mz, isNL
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
createAnnotatedSpectrum <- function(psm){

  ionAnnot <- unlist(strsplit(as.character(psm$Matches), ";"))
  ionType <- substr(ionAnnot,1,1)
  fragmentNb <- as.integer(gsub("[-(].*","",substr(ionAnnot,2,nchar(ionAnnot))))

  intensity <- as.numeric(unlist(strsplit(as.character(psm$Intensities), ";")))
  adjustedIntensity <-  intensity/psm$`Precursor Apex Fraction`
  mz <- as.numeric(unlist(strsplit(as.character(psm$Masses), ";")))
  mzDelta = as.numeric(unlist(strsplit(as.character(psm$'Mass Deviations [Da]'), ";")))

  charge <- rep(1,length(mz))
  isNL <- grepl("\\-",ionAnnot)
  isSelected <- rep(F,length(isNL))

  # heavy / light
  labelState <- grepl("Lys8|Arg10",psm$Modifications)

  ### non +1 fragment charge states
  higherCS <- grepl("\\+",ionAnnot)
  if(T %in% higherCS){
    chargeTmp <- gsub(".*\\(","",ionAnnot)
    chargeTmp <- gsub("\\+\\)","",chargeTmp)
    charge[higherCS] <- as.numeric(chargeTmp[higherCS])
  }

  # theoretical fragment mz
  mz = mz + mzDelta

  # ret <- list()
  # ret$spectrum <- data.frame(intensity,adjustedIntensity,ionAnnot,ionType,fragmentNb,charge,mz,isNL,isSelected)
  # ret$precMz <- (psm$Mass+psm$Charge*1.007316)/psm$Charge
  # ret$precIntensity <- psm$`Precursor Intensity`
  # ret$isHeavy <- labelState
  # ret$retentionTime <- psm$`Retention time`
  # ret$iRT <- suppressWarnings(ifelse(is.null(psm$iRT),NA,psm$iRT))
  # ret$charge <- as.numeric(psm$Charge)
  # ret$protein <- as.character(psm$Proteins)
  # ret$peptide <- as.character(psm$Sequence)
  # ret$mqResIdx <- rownames(psm)
  #
  # return(ret)

  spectrum = data.frame(intensity,adjustedIntensity,ionAnnot,ionType,fragmentNb,charge,mz,isNL,isSelected
                        , precMz = (psm$Mass+psm$Charge*1.007316)/psm$Charge
                        , precIntensity =  psm$`Precursor Intensity`
                        , isHeavy = labelState
                        , retentionTime = psm$`Retention time`
                        , iRT = suppressWarnings(ifelse(is.null(psm$iRT),NA,psm$iRT))
                        , precCharge = as.numeric(psm$Charge)
                        , protein = as.character(psm$Proteins)
                        , peptide = as.character(psm$Sequence)
                        , mqResIdx = rownames(psm)
                        )
  return(spectrum)
}

#' Get linear model predicting iRT as a function retention time (column name "rt")
#' @param tb tibble including column "Retention time" and "Sequence"
#' @return list including fit rlm object,  data.frame peptide,rt,irtRef
#' @export
#' @note  No note
#' @details Uses Robust Linear Regression. Retention times of at least 3 iRT unique peptides needs to be listed in the input tibble
#' @references Using iRT, a normalized retention time for more targeted measurement of peptides, Escher et al. (2012), \url{http://www.ncbi.nlm.nih.gov/pubmed/22577012}
#' @examples print("No examples")
getIRTModel <- function(tb){

  irtDF <- data.frame(peptide = tb$Sequence,rt=tb$`Retention time`,irtRef= IRTPEPTIDES[ match( tb$Sequence , rownames(IRTPEPTIDES)),] ) %>% na.omit

  ### require at least 3 identified irt peptides
  if(length(unique(irtDF$peptide)) < 3){
    stop("Not enough irt peptides detected")
  }

  ret = list()
  ret$fit = rlm(irtRef ~ rt ,data=irtDF)
  ret$irtPeptides = irtDF

  return(ret)

}

#' Add iRT metric to data.frame
#' @param tb tibble containing colum labelled "Retention time"
#' @param fit lm object
#' @return vector of normalised rt (empirical irt)
#' @export
#' @note  No note
#' @details No details
#' @references Using iRT, a normalized retention time for more targeted measurement of peptides, Escher et al. (2012), \url{http://www.ncbi.nlm.nih.gov/pubmed/22577012}
#' @examples print("No examples")
getEmpiricalIRT <- function(tb,fit){
  return(predict(fit, newdata=cbind(tb,rt = tb$`Retention time` )))
}

#' Plot IRT calibration curve
#' @param irtModel list
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotIRTCalibration<- function(irtModel
                              , cex=1.5
                              , col="blue"
                              , pch=19
                              ,...){

  fit = irtModel$fit
  irt = irtModel$irtPeptides
  irt <- irt[irt$peptide %in% rownames(IRTPEPTIDES),]

  plot(irt$rt , IRTPEPTIDES[match(as.character(irt$peptide),rownames(IRTPEPTIDES)),]
       , xlab="Retention Time"
       , ylab="iRT Reference"
       , cex.lab=cex
       , cex.axis=cex
       , cex=cex
       , col=col
       , pch=pch
       , ...)

  # missing irt
  abline(h=setdiff(IRTPEPTIDES$refIRT,IRTPEPTIDES[match(as.character(irt$peptide),rownames(IRTPEPTIDES)),] ), col ="red", lwd=1.5)

  rSquared <- round(cor(irt$rt , IRTPEPTIDES[as.character(irt$peptide),1],use="complete")^2,2)
  legend("topleft",as.expression(c(bquote(R^2*"" == .(rSquared)))), box.lwd=0,box.col=0,cex=cex)
  legend("bottomright", paste("iRT = ",round(coef(fit)[2],1)," x rt + ",round(coef(fit)[1],1),sep="") , box.lwd=0,box.col=0,cex=cex)
  abline(fit,lwd=1.5)
}

#' Get Fragment Sub Sequnece
#' @param peptide character string
#' @param ionType c("a","b","x","y")
#' @param fragmentNb  integer
#' @return character string
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getFragmentSequence <- function(peptide=peptide,ionType=ionType,fragmentNb=fragmentNb){

  # cast
  peptide = peptide %>% as.character

  if(ionType %in% c("y","x")){
    substr(peptide,(nchar(peptide)-fragmentNb)+1,nchar(peptide))
  }else if((ionType %in% c("a","b"))){
    substr(peptide,1,fragmentNb)
  }else{
    warning("Unknown ionType ", ionType)
    return(NA)
  }
}


#' Get Mz Shift of complementary peptide/fragment ion
#' @param aaSeq character string
#' @param charge default 1
#' @param isHeavy c(T,F)
#' @return numeric
#' @export
#' @note  No note
#' @details Calculate mass shift of complimentary spectrum.
#'  \itemize{
#'  \item PETIDEK (light) -> 8.014199
#'  \item PETIDEK (heavy) -> -8.014199
#'  \item PETIDER (light) -> 10.008269
#'  \item PETIDER (heavy) -> -10.008269
#'  }
#' @references NA
#' @examples print("No examples")
getLabelMzShift <- function(aaSeq, charge=1,isHeavy=T){

  mzShift = 0
  kShift <- 8.014199
  rShift <- 10.008269

  ### ONLY LABELLED AT THE C-TERMINI
  mzShift <- 0
  if(grepl("K$" ,aaSeq)){
    mzShift <- kShift / charge
  }else if(grepl("R$" ,aaSeq)){
    mzShift <- rShift / charge
  }

  if(isHeavy){
    return(-mzShift)
  }
  return(mzShift)
}

#' Created complimentary Heavy/Light annotated spectrum with updates precursor and fragment m/z values.
#' @param annotSpectrum Annotated Spectrum object
#' @return annotSpectrum Annotated Spectrum object
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getComplementaryLabelSpectrum <- function(annotSpectrum){

  # update precursor mz
  annotSpectrum$precMz <-  annotSpectrum$precMz + getLabelMzShift(aaSeq = annotSpectrum$peptide[1], charge=annotSpectrum$charge[1], isHeavy=annotSpectrum$isHeavy[1])

  # update fragment mz
  for(i in 1:nrow(annotSpectrum)){
    aaSeq <- getFragmentSequence(peptide=annotSpectrum$peptide[1]  ,ionType=as.character(annotSpectrum[i,]$ionType),fragmentNb=as.numeric(annotSpectrum[i,]$fragmentNb))
    annotSpectrum[i,]$mz = annotSpectrum[i,]$mz + getLabelMzShift(aaSeq = aaSeq, charge=as.numeric(annotSpectrum[i,]$charge), isHeavy=annotSpectrum$isHeavy[1])

  }
  annotSpectrum$isHeavy <- !annotSpectrum$isHeavy

  return(annotSpectrum)
}

#' Digest protein
#' @param proteinSeq protein sequence
#' @param proteaseRegExp protease Regular Expression
#' @param nbMiscleavages default 0
#' @return vector of peptides
#' @export
#' @note  No note
#' @details No details
#' @examples print("No examples")
getPeptides <- function(proteinSeq,proteaseRegExp="[KR](?!P)",nbMiscleavages=0){

  allAA <- as.vector(unlist(strsplit(proteinSeq,"")))

  ### get all full cleaved peptides without cleaved residue
  fcPeptides <- strsplit(proteinSeq,proteaseRegExp,perl=TRUE)[[1]]

  ### add cleaved residue
  matchPos <- gregexpr(proteaseRegExp,proteinSeq,perl=TRUE)[[1]]
  separator <- allAA[ matchPos ]
  ###if c-term peptide isn't tryptic
  if(length(separator) < length(fcPeptides)) separator <- c(separator,"")
  fcPeptides <- paste(fcPeptides,separator,sep="")
  fcPeptides <- fcPeptides[nchar(fcPeptides) > 0 ]

  ### if no mis-cleavages that's it
  if(nbMiscleavages == 0){
    return(fcPeptides)
  }

  allPeptides <- c()
  ### handle miscleavages
  for(i in 1:length(fcPeptides)){
    #cat(fcPeptides[i]," ",i," -------------\n")
    for(j in i:(i+nbMiscleavages)){
      if(j <= length(fcPeptides)){
        pept <- paste(fcPeptides[i:j],collapse="")
        #cat(j," ---",pept,"\n")
        allPeptides <- c(allPeptides,pept)
      }
    }
  }

  return(allPeptides[nchar(allPeptides) > 0])

}

#' Get peptide candidates per protein, appying some hardcoded files
#' @param proteins list
#' @param dispProgressBar T/F
#' @param peptideLengthRange integer vector default c(6,21)
#' @param proteaseRegExp character '[KR](?!P)' - trypsin
#' @param nbMiscleavages 0
#' @param trimAC TRUE/FALSE default TRUE sp|P62258|1433E_HUMAN -> P62258
#' @param exclusivePeptides TRUE/FALSE get exclusive peptides only i.e. peptides mapping to a single protein. default T
#' @return data.frame peptide,protein
#' @export
#' @note  No note
#' @references NA
#' @import Biobase
#' @examples print("No examples")
digestProteome = function(proteins, proteaseRegExp="[KR](?!P)",dispProgressBar=T, peptideLengthRange=c(6,21), nbMiscleavages=0, exclusivePeptides=T, trimAC=T  ){

  # progress bar
  if(dispProgressBar) cat("Digesting Proteome \n")
  pbSum <- txtProgressBar(min = 0, max = length(proteins), style = 3)

  peptideDf = data.frame()

  for(i in 1:length(proteins)){

    prot <- names(proteins)[i]
    ### increment progressbar
    if(dispProgressBar) setTxtProgressBar(pbSum, i)
    peptides <- getPeptides(unlist(proteins[prot]), proteaseRegExp= proteaseRegExp, nbMiscleavages=nbMiscleavages)

    peptides = peptides[nchar(peptides) %in%  peptideLengthRange[1]:peptideLengthRange[2] ]

    if(length(peptides) > 0 ){
      peptideDf = rbind(peptideDf,  cbind(peptides,  prot) )
    }
  }

  # close progress bar
  setTxtProgressBar(pbSum, i)
  close(pbSum)

  # rename columns
  if(nrow(peptideDf) > 0)  names(peptideDf) = c("peptide","protein")

  # keep exclusive peptides only (i.e. peptides matching only one protein)
  if(exclusivePeptides){peptideDf = group_by(peptideDf, peptide) %>% filter(length(protein) == 1) %>% data.frame }

  # sp|P62258|1433E_HUMAN -> P62258
  # sp|Q9JII5|DAZP1_MOUSE -> Q9JII5
  if(trimAC) peptideDf$protein = gsub("(^[a-z]{2}\\|)([A-Z0-9\\-]*)(\\|.*)","\\2",peptideDf$protein)

  peptideDf$length = peptideDf$peptide %>% as.character %>% nchar

  return(peptideDf)

}

