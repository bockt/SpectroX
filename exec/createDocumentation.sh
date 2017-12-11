#!/bin/sh
Rscript /Users/ahrnee-adm/dev/R/workspace/SpectroX/exec/roxygenize.R
R CMD Rd2pdf /Users/ahrnee-adm/dev/R/workspace/SpectroX --output=/Users/ahrnee-adm/dev/R/workspace/SpectroX//inst/manuals/SpectroX-man.pdf --force --no-preview
