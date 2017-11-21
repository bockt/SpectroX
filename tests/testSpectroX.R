

library(tidyverse)





#' Parse MaxQuant msms.txt
#' @param file path
#' @param pepThrs numeric default 0.05
#' @param targetPeptides character default NA
#' @param targetProteins character default NA
#' @param filterContaminants TRUE
#' @param contaminantRegExp '^CON_'
#' @param selectedPTMRegExp NA
#' @param filterNonExclusivePeptides default TRUE
#' @return data.frame of maxQuant psm level results
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
parseMqxQuantMSMS = function(file
                             , pepThrs=0.05
                             , targetPeptides=NA
                             , targetProteins=NA
                             , filterContaminants = T
                             , contaminantRegExp = '^CON_'
                             , selectedPTMRegExp = NA
                             , filterNonExclusivePeptides = T
                             ){
  # do not filter out non exclusive peptides if target peptide list is provided
  if(length(targetPeptides) > 0) filterNonExclusivePeptides = F

  # do not filter out iRT peptides when filtering put contaminants

  # discard decoy protein matches

  # Type filter
  # The type of the feature.
  # 'MSMS' – for an MS/MS spectrum without an MS1 isotope pattern assigned.
  # 'ISO-MSMS' – MS1 isotope cluster identified by MS/MS.
  # 'MULTI-MSMS' – MS1 labeling cluster identified by MS/MS.
  # 'MULTI-SECPEP' – MS1 labeling cluster identified by MS/MS as second peptide.
  # 'MULTI-MATCH' – MS1 labeling cluster identified by matching between runs.
  # In case of label-free data there is no difference between 'MULTI' and 'ISO'.

}


#file = "inst/msms.txt"
#file = "/Users/ahrnee-adm/dev/R/workspace/SpectroUtils/inst/testData/msms.tsv"
file = "/Users/ahrnee-adm/dev/R/workspace/SpectroUtils/inst/testData/old/msmsTest.tsv"


tb = read_tsv(file)

dim(tb)



table(tb$id) %>% max


boxplot(tb$`Delta score` ~ tb$Type)

class(tb)

str(tb)
names(tb)

table(tb$id)

nrow(tb)

table(tb$Modifications)
table(tb$Type)

tb$`Isotope index`

tb$Type




subset(tb, Type == "MSMS")$"Precursor Apex Fraction" %>% max

subset(tb, Type == "MULTI-SECPEP")$"Precursor Apex Fraction" %>% max

subset(tb, Type == "MULTI-MSMS")$"Precursor Apex Fraction" %>% max

tb$"Precursor Apex Fraction" %>% hist


tb2 = read_tsv("../SpectroUtils/inst/testData/msms.tsv")

setdiff(names(tb),names(tb2)  )
setdiff(names(tb2),names(tb)  )

tb2$Precursor

