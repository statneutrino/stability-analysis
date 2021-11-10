args = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

print(args)

if (is.null(args$input)){
    input = './run_tms_steroid_stability.Rmd'
} else {
    input = args$input
}
if (is.null(args$filename)){
    filename = "tms_steroid_fitting.html"
} else {
    filename = args$filename
}
if (is.null(args$dir)){
    dir = "./results"
} else {
    dir = args$dir
}

Sys.setenv(RSTUDIO_PANDOC='C:/Program Files/RStudio/bin/pandoc')
rmarkdown::render(
    input, 
    output_file = filename,
    output_dir = dir)