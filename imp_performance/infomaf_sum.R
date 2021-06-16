#!/usr/bin/env Rscript

######################
## Annotate Results ##
######################

#Author: Andres Kaufmann


library(data.table)

options(stringsAsFactors = F)
args <- commandArgs(TRUE)

#Print help as default

if(length(args) < 1) {
  args <- c("--help")
}

print.help <- function(){
  cat("
      Wellcome
      
      This R script binds the results from imputation in a single file and annotates the ancestry of each
      position according to local ancestry results.
      
      Arguments:
      
      
      --dir= The directory of the results
      --ind= The index of the individual (number)
      --out = The output directory
      
      
      Good luck with you research.
      ")
  
  q(save = "no")
}

if("--help" %in% args){
  print.help()
}

###############################################
### Parse arguments in the form --arg=value ###
###############################################

parse_args <- function(x) strsplit(sub("^--", "", x), "=")
argsdf <- as.data.frame(do.call("rbind", parse_args(args)))
argsL <- as.list(as.character(argsdf$V2))
names(argsL) <- argsdf$V1

##########################################
### Check default arguments are parsed ###
##########################################

if(is.null(argsL$ind)){
  print.help()
  stop("Your individual")
} else {
  ind <- argsL$ind
}
if(is.null(argsL$dir)){
  print.help()
  stop("You have to provide a directory with the results")
} else{
  dir <- argsL$dir
}
if(is.null(argsL$out)){
  print.help()
  stop("You have to provide an output directory")
} else {
  out <- argsL$out
}

ind <- as.numeric(ind)

####################
### START SCRIPT ###
####################
#Defining directories and variables
impute_dir <- paste0(dir, "/", "indiv", ind)
amr_ids <- fread("/well/hill/adriank/MEX_BImp/data/array_positions/AMR_ind.txt")


#Define the AMR id and read the .bed files
id <- amr_ids[ind]$ID
pop <- amr_ids[ind]$POP

info_files <- list.files(path = impute_dir, pattern = "_info$", full.names = T)
info_dat <- lapply(info_files, fread)
info_dat <- rbindlist(info_dat)
#Crear MAF
info_dat[, MAF := ifelse(exp_freq_a1 <= 0.5, exp_freq_a1, 1-exp_freq_a1)]
inf_filt <- nrow(info_dat[info >= 0.3])
maf_filt <- nrow(info_dat[MAF >= 0.01])
two_filt <- nrow(info_dat[info >= 0.3 & MAF >=0.01])
filters <- data.table(INFO = inf_filt, MAF = maf_filt, ALL = two_filt)
subset <- info_dat[, c("info", "MAF")]

#################
####Write data###
#################
outdir <- paste0(out, "/", pop, "/")
fwrite(filters, file = paste0(outdir, "ind",ind, "_filters.txt"), sep = "\t")
fwrite(subset, file = paste0(outdir,"ind",ind, "_infomaf.txt"), sep = "\t")