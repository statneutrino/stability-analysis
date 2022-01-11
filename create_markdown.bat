@echo off
PATH "C:\Program Files\R\R-4.0.5\bin\x64";%path%
@echo on
Rscript --vanilla create_report.R --input="produce_tables_urine.Rmd" --filename="urine_output_tables.html" --dir="results"