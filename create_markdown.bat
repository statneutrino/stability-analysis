@echo off
PATH "C:\Program Files\R\R-4.0.5\bin\x64";%path%
@echo on
Rscript -e "Sys.setenv(RSTUDIO_PANDOC='C:/Program Files/RStudio/bin/pandoc'); rmarkdown::render('./run_tms_steroid_stability.Rmd')"