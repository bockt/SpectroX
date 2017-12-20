suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
VERSION = 0.0

source( "/Users/ahrnee-adm/dev/R/workspace/SpectroX/R/UserOptions.R")

### USER CMD LINE OPTIONS
cmdlineOptions <- getCMDLineOptions(version=VERSION)
uO = getUserOptions(cmdlineOptions)
### USER CMD LINE OPTIONS END




