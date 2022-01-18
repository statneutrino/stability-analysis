@echo off
PATH "C:\Program Files\R\R-4.1.2\bin";%path%
@echo on
Rscript --vanilla create_report.R --input="produce_tables_urine.Rmd" --filename="output_tables_urine.html" --dir="results"