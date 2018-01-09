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
#' @param pepCutoff numeric default 0.05
#' @param targetPeptides character default NA
#' @param targetProteins character default NA
#' @param filterContaminants TRUE
#' @param contaminantRegExp '^CON_'
#' @param selectedPTMRegExp NA
#' @param filterNonExclusivePeptides default TRUE
#' @param pepLength default [1,Inf)
#' @param chargeState default [1,Inf)
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
                             , pepCutoff=0.05
                             , targetPeptides=NA
                             , targetProteins=NA
                             , filterContaminants = T
                             , contaminantRegExp = '^CON_'
                             , selectedPTMRegExp = NA
                             , filterNonExclusivePeptides = T
                             , pepLength = c(1,Inf)
                             , chargeState = c(1,Inf)
                             , label = NA
                             , maxMissedCleavages = 0
                             , keepBestSpectrumOnly = T
                             , requiredPepSeqRegExp = "."
                             ,...

){

  tb = suppressMessages(read_tsv(file))

  # is irt peptide (unmodified)
  tb$isIRT = paste0("_",rownames(IRTPEPTIDES),"_" , collapse="|") %>% grepl(.,tb$`Modified sequence`)

  # rename irt protein entry
  #tb$Proteins[tb$isIRT] %>% print
  tb$Proteins[tb$isIRT] = "BiognosysIRT"

  # target peptide filter
  if(!is.na(targetPeptides[1])){
    # targetPeptides = c("AAAAAAAAALGVR", "AAAAADEDEEADEEDEAEAGGMGPQAR" )
    tb = subset(tb, (Sequence %in% targetPeptides) | isIRT  )
    # disable additional peptide filters
    filterNonExclusivePeptides = F
    pepLength = c(1,Inf)
    maxMissedCleavages = Inf
    requiredPepSeqRegExp = "."
  }
  #Filter out non exclusive peptides
  if(filterNonExclusivePeptides){
    tb = subset(tb, !grepl("\\;",Proteins))
  }else if(!is.na(targetPeptides[1])){ # warn if target peptides are not exclusive
    nonExclPeptides =  subset(tb, grepl("\\;",Proteins))
    if(nrow(nonExclPeptides) > 0 ){
      cat("Warning - Found Non-Exclusive Peptides:")
      nonExclPeptides[c("Sequence","Proteins")] %>% as.data.frame() %>% unique %>% print
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
  keep = keep & (tb$PEP <  pepCutoff)
  # decoy filter
  keep = keep & (is.na(tb$Reverse) | tb$Reverse != "+" )
  # peptide length
  keep = keep &  ((tb$Length >= min(pepLength)) & (tb$Length <= max(pepLength)) )


  # keep charge state
  keep = keep &  ((tb$Charge >= min(chargeState)) & (tb$Charge <= max(chargeState)) )
  # keep modifs
  #selectedPTMRegExp = "Oxidation"
  # if(!is.na(selectedPTMRegExp)){
  #   keep =  keep & (paste(c("Unmodified","Arg10","Lys8",selectedPTMRegExp),collapse="|") %>% grepl(.,tb$Modifications))
  # }else{
  #   keep =  keep & (paste(c("Unmodified","Arg10","Lys8"),collapse="|") %>% grepl(.,tb$Modifications))
  # }

  # get all modifications seached for
  excludedPTM = getSearchedModifications(tb,...)

  cat("INFO: allowed PTMs:", ifelse(is.na(excludedPTM[grepl(selectedPTMRegExp,excludedPTM)][1]),"None", paste(excludedPTM[grepl(selectedPTMRegExp,excludedPTM)],collapse=":")),"\n")
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
  #keep = keep | ((paste(IRTPEPTIDES, collapse="|") %>% grepl(.,tb$`Modified sequence`)) &  (tb$PEP <  pepCutoff))
  keep = keep | (tb$isIRT & (tb$PEP <  pepCutoff))

  # apply filter
  tb = subset(tb, keep)

  # filterBestSpectrum
  # keep top-scoring spectrum only, per peptidoform
  if(keepBestSpectrumOnly) tb = group_by(tb, Sequence, Modifications) %>% filter(rank(Score,ties.method = "first") == length(Score))

  # casting
  suppressWarnings(tb$`Precursor Apex Fraction` %<>% as.numeric)
  suppressWarnings(tb$`Precursor Intensity` %<>% as.numeric)

  return(tb)

}

#' Get list of variable modifiections condsidered by MaxQuant
#' @param df tible or data.frame
#' @param ignoreArgLysIsoLabel default T Do not consider Arg Lys C,N heavy isotope label
#' @return character vector of modification names (does not return "Arg10","Lys8" labels )
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getSearchedModifications = function(tb, ignoreArgLysIsoLabel=T){

  mods = names(tb)[grepl("Probabilities$",names(tb))] %>% gsub(" Probabilities$","",.)
  # add n-term mods
  mods = c(mods,names(tb)[grepl("N\\-term",names(tb))])
  # remove "Arg10"         "Lys8"
  if(ignoreArgLysIsoLabel) mods = mods[!(mods %in% c("Arg10","Lys8"))]
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
createLibrarySpectrum <- function(psm){

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
                        , precApexIntensity = psm$`Precursor Intensity` / psm$`Precursor Apex Fraction`
                        , psmScore = psm$Score
                        , protein = as.character(psm$Proteins)
                        , proteinDescription =psm$`Protein Names`
                        , peptide = as.character(psm$Sequence)
                        , ptm = psm$Modifications
                        , mqResIdx = psm$`Scan number`
                        , isIRT = psm$isIRT


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
# getLabelMzShift <- function(aaSeq, charge=1,isHeavy=T){
#
#   mzShift = 0
#   kShift <- 8.014199
#   rShift <- 10.008269
#
#   ### ONLY LABELLED AT THE C-TERMINI
#   mzShift <- 0
#   if(grepl("K$" ,aaSeq)){
#     mzShift <- kShift / charge
#   }else if(grepl("R$" ,aaSeq)){
#     mzShift <- rShift / charge
#   }
#
#   if(isHeavy){
#     return(-mzShift)
#   }
#   return(mzShift)
# }

#' Created complimentary Heavy/Light annotated spectrum with updates precursor and fragment m/z values.
#' @param annotSpectrum Annotated Spectrum object
#' @return annotSpectrum Annotated Spectrum object
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
# getComplementaryLabelSpectrum <- function(annotSpectrum){
#
#   # update precursor mz
#   annotSpectrum$precMz <-  annotSpectrum$precMz + getLabelMzShift(aaSeq = annotSpectrum$peptide[1], charge=annotSpectrum$charge[1], isHeavy=annotSpectrum$isHeavy[1])
#
#   # update fragment mz
#   for(i in 1:nrow(annotSpectrum)){
#     aaSeq <- getFragmentSequence(peptide=annotSpectrum$peptide[1]  ,ionType=as.character(annotSpectrum[i,]$ionType),fragmentNb=as.numeric(annotSpectrum[i,]$fragmentNb))
#     annotSpectrum[i,]$mz = annotSpectrum[i,]$mz + getLabelMzShift(aaSeq = aaSeq, charge=as.numeric(annotSpectrum[i,]$charge), isHeavy=annotSpectrum$isHeavy[1])
#
#   }
#   annotSpectrum$isHeavy <- !annotSpectrum$isHeavy
#
#   return(annotSpectrum)
# }

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

#' Get peptide candidates per protein
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
#' @examples print("No examples")
digestProteome = function(proteins, proteaseRegExp="[KR](?!P)",dispProgressBar=T, peptideLengthRange=c(6,21), nbMiscleavages=0, exclusivePeptides=T, trimAC=T  ){

  # progress bar
  if(dispProgressBar) cat("Digesting Proteome \n")
  pbSum <- txtProgressBar(min = 0, max = length(proteins), style = 3)

  nbProteins =length(proteins)
  # assume 100 peptides per protein
  allPeptides = rep(NA,nbProteins*100)
  allProteins =  rep(NA,nbProteins*100)

  j = 1
  for(i in 1:nbProteins){

    prot <- names(proteins)[i]
    ### increment progressbar
    if(dispProgressBar) setTxtProgressBar(pbSum, i)

    # digest
    peptides <- getPeptides(unlist(proteins[prot]), proteaseRegExp= proteaseRegExp, nbMiscleavages=nbMiscleavages)
    # keep only peptides in accetped length range

    # add peptides and protein to vector
    k = (j+length(peptides)-1)
    allPeptides[j:k] = peptides
    allProteins[j:k] = prot
    j = k+1

  }

  # create df
  peptideDf = na.omit(data.frame(allPeptides,allProteins))
  names(peptideDf) = c("peptides","proteins")

  # filter by peptide length
  pepLength = nchar(peptideDf$peptides %>% as.character)
  peptideDf = subset(peptideDf,  (pepLength >=  min(peptideLengthRange)) & (pepLength <= max(peptideLengthRange)))

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


#'# Export top X peptide per protein (ordered by 'adjustedIntensitySum') and add theoretical peptides ranked by length if fewer than X peptides were identified by MaxQuant
#' @param spectralLibrary tibble
#' @param theoPeptides data.frame 'peptide' 'protein' 'length'
#' @param targetProteins list of protein accesion numbers e.g. P62258. default all proteins in spectralLibrary
#' @param nbPeptidesPerProtein number of peptides per protein to be exported
#' @param outFile path to output file defult tmp.xls in tmp dir
#' @export
#' @note  No note
#' @references NA
#' @examples print("No examples")
proteotypicPeptideExport = function(spectralLibrary
                                    ,theoPeptides=NA
                                    ,targetProteins=unique(spectralLibrary$protein)
                                    , nbPeptidesPerProtein=5
                                    , outFile= paste0(tempdir(),"/tmp.xls")){

  # output protein, peptide, rankingMetric, rank (if multiple charge states pick most intense)
  xlsOut = spectralLibrary[c("protein","peptide","ptm","rankingMetric","isIRT","adjustedIntensitySum", "psmScore", "precApexIntensity")]%>%
    unique %>% group_by(protein,peptide,ptm) %>% top_n(1,rankingMetric) %>%
    split(.$protein) %>%
    map_df(~ data.frame(peptide=.x$peptide
                        , ptm = .x$ptm
                        , protein = .x$protein
                        , isIRT = .x$isIRT
                        , adjustedIntensitySum =  .x$adjustedIntensitySum
                        , psmScore =  .x$psmScore
                        , precApexIntensity =  .x$precApexIntensity
                        , rankingMetric= .x$rankingMetric
                        , rank = length(.x$rankingMetric) - rank(.x$rankingMetric)+1
    ))

  # get proteins without enough peptides
  # nbMissingPeptidesPerProt = subset(xlsOut,rank <= nbPeptidesPerProtein) %>%
  #   group_by(. ,protein) %>%

  nbMissingPeptidesPerProt = group_by(xlsOut ,protein) %>%
  summarise(nbMissingPep = nbPeptidesPerProtein - max(rank,na.rm=T)) %>%
  subset(nbMissingPep > 0) %>% as.data.frame
  rownames(nbMissingPeptidesPerProt) = nbMissingPeptidesPerProt$protein

  # add proteins with no identified peptides
  missingProteins =  setdiff(targetProteins, nbMissingPeptidesPerProt$protein)
  if(sum(!is.na(missingProteins)) > 0 ){
    nbMissingPeptidesPerProt =rbind(nbMissingPeptidesPerProt
                                    ,data.frame(protein=missingProteins,nbMissingPep =nbPeptidesPerProtein, row.names=missingProteins)
    )
  }

  # add theo peptides
  if(!is.null(theoPeptides %>% dim)){
    for(prot in nbMissingPeptidesPerProt$protein ){
      # if protein has theo peptides
      if(prot %in% theoPeptides$protein ){
        subTP = subset(theoPeptides, protein %in% prot )
        # order peptides by length
        subTP = subTP[ order(subTP$length),]
        # add nbMissingPep peptides
        sel = 1:min(nbMissingPeptidesPerProt[match(prot, rownames(nbMissingPeptidesPerProt)),]$nbMissingPep, nrow(subTP))
        xlsOut = rbind(xlsOut, cbind(subTP[sel,c("peptide", "protein")]
                                     , rankingMetric=NA
                                     , ptm =NA
                                     , isIRT = F
                                     , rank=nbPeptidesPerProtein
                                     , adjustedIntensitySum =  NA
                                     , psmScore =  NA
                                     , precApexIntensity =  NA

                                     ) )
      }
    }
  }


  #keep all irts
  write.table(file=outFile,  subset(xlsOut, (rank <= nbPeptidesPerProtein) | isIRT ) %>% arrange(.,protein,rank) ,row.names=F,sep = "\t")
  cat("CREATED proteotypicPeptideExport FILE: ", outFile,"\n")

}




#' Create Spectral Library
#' @param tb tibble maxQuant spectrum level search results
#' @param minFragNb min frag number default 3 (i.e. b3 and y3 will be kept)
#' @param minNbTransitions minimum number of fragments
#' @param minBasePeakFraction minimum intnsity fraction of base peak (most intense kept fragment)
#' @param ionTypeFilter selected ion type default a,b,x,y
#' @param includeNL include neutral loss peaks defualt TRUE
#' @param rankingMetric character. column name, used for ranking peptide, default "adjustedIntensitySum", c("adjustedIntensitySum","psmScore","precApexIntensity")
#' @return data.frame
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
createSpectralLibrary = function(tb, minFragNb = 3, minNbTransitions = 5,maxNbTransitions = 5, minBasePeakFraction = 0, ionTypeFilter = c("a","b","x","y"), includeNL = T, rankingMetric = "adjustedIntensitySum"  ){

  cat("Parsing Spectra\n")
  pbSum <- txtProgressBar(min = 0, max = nrow(tb), style = 3)
  spectralLibrary = data.frame()
  for(rownb in 1:nrow(tb) ){

    setTxtProgressBar(pbSum, rownb)
    annotSpec = createLibrarySpectrum(tb[rownb,])

    ### apply spectrum filters
    # fragment number filter
    keep = with(annotSpec,  (fragmentNb >= minFragNb) & grepl(paste(ionTypeFilter,collapse="|"), annotSpec$ionType) )
    # neutral loss peak filter
    if(!includeNL) keep = keep &  !annotSpec$isNL
    annotSpec =   subset(annotSpec,keep )

    #fragment intensity as fraction of base peak
    annotSpec$basePeakIntFrac = annotSpec$intensity / suppressWarnings(max(annotSpec$intensity,na.rm=T))
    # filter minimum base peak intensity fraction
    annotSpec = subset(annotSpec,basePeakIntFrac >= minBasePeakFraction )

    # make sure enough fragments
    if(nrow(annotSpec) >= minNbTransitions){

      #select maxNbTransitions fragments
      annotSpec = arrange(annotSpec, intensity)[1:min(nrow(annotSpec),maxNbTransitions),]

      # calculate relative intensity based on fragment seleciton
      annotSpec$relativeIntensity = (annotSpec$intensity / max(annotSpec$intensity))*100

      # sum of adjusted intensity (basis for ranking)
      annotSpec$adjustedIntensitySum =  annotSpec$adjustedIntensity %>% sum(.,na.rm=T)

      # default ranking metric - summed adjusted intensity
      if(rankingMetric %in% names(annotSpec) & is.numeric(annotSpec[,rankingMetric])){
        annotSpec$rankingMetric =  annotSpec[,rankingMetric][1]
      }else{
        stop("createSpectralLibrary: Invalid rankingMetric", rankingMetric)
      }

      #add to library, SLOW append
      spectralLibrary = rbind(spectralLibrary,annotSpec )
    }
  }
  # close progress bar
  setTxtProgressBar(pbSum, rownb)
  close(pbSum)

  return(spectralLibrary)

}

#' Create complementary isotope (Arg, Lys H/L) spectral library
#' @param sl tibble or data frame
#' @return data.frame
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
createComplementaryIsotopeLibrary = function(sl){

  kShift =8.014199
  rShift = 10.008269

  hasKCTerm =  grepl("K$" ,sl$peptide)
  hasRCTerm =  grepl("R$" ,sl$peptide)

  # correct precursor m/z
  # heavy Lys c-term peptide
  sl$precMz[sl$isHeavy] = sl$precMz[sl$isHeavy] - (ifelse(hasKCTerm,kShift,0)[sl$isHeavy] / sl$precCharge[sl$isHeavy] )
  # heavy Arg c-term peptide
  sl$precMz[sl$isHeavy] = sl$precMz[sl$isHeavy] - (ifelse(hasRCTerm,rShift,0)[sl$isHeavy] / sl$precCharge[sl$isHeavy] )
  # light Lys c-term peptide
  sl$precMz[!sl$isHeavy] = sl$precMz[!sl$isHeavy] + (ifelse(hasKCTerm,kShift,0)[!sl$isHeavy] / sl$precCharge[!sl$isHeavy] )
  # light Arg c-term peptide
  sl$precMz[!sl$isHeavy] = sl$precMz[!sl$isHeavy] + (ifelse(hasRCTerm,rShift,0)[!sl$isHeavy] / sl$precCharge[!sl$isHeavy] )

  # correct c-term fragment m/z
  isCTermFrag =  sl$ionType %in% c("x","y","z")
  # heavy Lys c-term fragment
  sl$mz[sl$isHeavy & isCTermFrag] = sl$mz[sl$isHeavy & isCTermFrag] - (ifelse(hasKCTerm,kShift,0)[sl$isHeavy & isCTermFrag] / sl$charge[sl$isHeavy & isCTermFrag] )
  # heavy Arg c-term fragment
  sl$mz[sl$isHeavy & isCTermFrag] = sl$mz[sl$isHeavy & isCTermFrag] - (ifelse(hasRCTerm,rShift,0)[sl$isHeavy & isCTermFrag] / sl$charge[sl$isHeavy & isCTermFrag] )
  # light Lys c-term fragment
  sl$mz[!sl$isHeavy & isCTermFrag] = sl$mz[!sl$isHeavy & isCTermFrag] + (ifelse(hasKCTerm,kShift,0)[!sl$isHeavy & isCTermFrag] / sl$charge[!sl$isHeavy & isCTermFrag] )
  # light Arg c-term fragment
  sl$mz[!sl$isHeavy & isCTermFrag] = sl$mz[!sl$isHeavy & isCTermFrag] + (ifelse(hasRCTerm,rShift,0)[!sl$isHeavy & isCTermFrag] / sl$charge[!sl$isHeavy & isCTermFrag] )



  sl$isHeavy = !sl$isHeavy

  return(sl)

}

#' Write SpectroDive compatible xls file
#' @param sl tibble or data frame
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
spectroDiveExport =  function(sl,file){
  out = data.frame(eg_proteinId = sl$protein
                   ,eg_proteinDescription = sl$proteinDescription
                   ,eg_assayCollectionType = "external"
                   ,eg_iRT = sl$iRT
                   ,t_relativeFragmentIonIntensity = sl$relativeIntensity
                   ,tg_aminoAcidSequenceWithASCIIMod = paste(sl$peptide,sl$ptm,sep="_")
                   ,eg_aminoAcidSequence = sl$peptide
                   ,tg_precursorCharge = sl$precCharge
                   ,t_fragmentIonCharge = sl$charge
                   ,t_fragmentIonSeries = sl$ionType
                   ,t_fragmentIonNumber = sl$fragmentNb
                   ,tg_Q1 = sl$precMz
                   ,t_Q3 = sl$mz

  )
  write.table(out,file=file,sep="\t",row.names=F)
  cat("CREATED spectroDiveExport  FILE: ", file,"\n")
}

#' Write Spectronaut compatible xls file
#' @param sl tibble or data frame
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
spectronautExport =  function(sl,file){
  out = data.frame(protein_id = sl$protein
                   ,Workflow = "SPIKE_IN"
                   ,iRT = sl$iRT
                   ,RelativeFragmentIntensity = sl$relativeIntensity
                   ,IntModifiedSequence = paste(sl$peptide,sl$ptm,sep="_")
                   ,StrippedSequence = sl$peptide
                   ,IsotopicLabel = ifelse(sl$isHeavy,"heavy","light")
                   ,PrecursorCharge = sl$precCharge
                   ,FragmentCharge = sl$charge
                   ,FragmentType = sl$ionType
                   ,FragmentNumber = sl$fragmentNb
                   ,PrecursorMz = sl$precMz
                   ,FragmentMz = sl$mz
                   ,ExcludeFromQuantification = ""

  )
  write.table(out,file=file,sep="\t",row.names=F)
  cat("CREATED spectronautExport  FILE: ", file,"\n")
}

#' Write skyline compatible xls file
#' @param sl tibble or data frame
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
skylineExport =  function(sl,file){
  out = data.frame(sl$precMz,
                   sl$mz,
                   25,
                   sl$peptide,
                   sl$protein,
                   paste0(sl$ionType,sl$fragmentNb),
                   "",
                   ifelse(sl$isHeavy,"heavy","light")
  )
  write.table(out,file=file,sep="\t",row.names=F, col.names = F)
  cat("CREATED skylineExport FILE: ", file,"\n")
}

#' Parse list of target peptides/proteins
#' @param file path
#' @return list() peptides, proteins
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
parseTargetsFile = function(file){
  tb = suppressMessages(read_csv(file))
  names(tb) = names(tb) %>% tolower

  proteins = tb[which(grepl("^protein[s]{0,1}$",names(tb)))]
  peptides = tb[which(grepl("^peptide[s]{0,1}$",names(tb)))]

  ret = list()
  ret$proteins = NA
  ret$peptides = NA

  # give priority to peptides
  if(length(peptides) > 0 ){
    ret$peptides = peptides %>% unlist%>% as.vector() %>% unique()
  }else if(length(proteins) > 0){
    ret$proteins = proteins %>% unlist%>% as.vector()  %>% unique()
  }else{
    stop("parseTargetsFile Invalid File format:",  file)
  }
  return(ret)
}

#' plot peptide count per protein
#' @param spectralLibrary data.frame
#' @param acLenTrunc integer default 12 "SOMEVERYLONGAC" -> "SOMEVERY.."
#' @param protLabCex default 0.9
#' @export
#' @note  No note
#' @references NA
#' @examples print("No examples")
barplotPetideCountPerProtein = function(spectralLibrary, protLabCex=0.9, acLenTrunc=12, col = "blue" , ... ){
  pepPerProt = (spectralLibrary[c("protein","peptide","ptm") ] %>% unique)[,"protein"] %>%table %>% sort(.,decreasing = T)
  protLab =names(pepPerProt)

  if(length(protLab) < 30){
    bp =barplot(pepPerProt,las=2, xaxt="n", ylab = "Unique Peptide Count", col =col,...)

    axis(1, labels = FALSE,tick=F)
    # truncate labels "SOMEVERYLONGAC" -> "SOMEVERY.."
    sel = nchar( protLab) > acLenTrunc
    protLab[sel] =  paste0(substr( protLab,1,acLenTrunc),"..")[sel]
    ### 45 degree labels
    text(bp, par("usr")[3], srt=45, adj = 1.1,
         labels =protLab, xpd = TRUE ,cex=protLabCex)
  }else{
    #axis(1, labels = TRUE,tick=T, cex.axis=cex.axis)
    #mtext("Protein",side =1, line=3, cex=cex.lab)
    h = hist(pepPerProt,breaks =  1:max(pepPerProt)
         #,cex.axis=cex.axis
         #,cex.lab=cex.lab
         ,xlab = "Unique Peptide Count Per Protein"
         ,col=col
         ,main=""
         ,...
         )
    axis(1, labels = TRUE,tick=T,...)

  }
}

#' plot ranking metric vs peptide (per protein) barplot
#' @param df data.frame
#' @param pepLenTrunc integer AFADAMEVIPSTLAENAGLNPISTVTELR -> AFADAMEVIPSTLAE..
#' @param pepLabCex default 0.7
#' @param rankingMetric character
#' @export
#' @note  No note
#' @references NA
#' @examples print("No examples")
barplotPeptidesPerProtein = function(df, pepLenTrunc=12,pepLabCex=0.7,rankingMetric="rankingMetric",...){

  # display modified peptide in red
  col = c("red","blue")[((df$ptm == "Unmodified") | (df$ptm == "Lys8") | (df$ptm == "Arg10")) +1]

  # avoid log(0) issue
  df$rankingMetric = df$rankingMetric +1

  # order by intensity decreasing L-R
  df = arrange(df,-rankingMetric)
  bp =  barplot(df$rankingMetric %>% log10
                , main = df$protein[1]
                , xaxt = "n"
                #, ylab="Adj. Fragment Int. Sum\n log10"
                ,ylab = paste0(rankingMetric,"\n log10")
                ,col=col
                ,...
  )

  # truncate labels AFADAMEVIPSTLAENAGLNPISTVTELR -> AFADAMEVIPSTLAE..
  df$peptide %<>% as.character
  sel = nchar( df$peptide) > pepLenTrunc
  df$peptide[sel] = paste0(substr( df$peptide,1,pepLenTrunc), sep="..")[sel ]
  df$peptide = paste0(df$peptide,"/",df$precCharge)
  ### 35 degree labels
  axis(1, labels = FALSE,tick=F)
  text(bp , par("usr")[3], srt=35, adj = 1.1,
       labels =df$peptide, xpd = TRUE ,cex=pepLabCex)

}
