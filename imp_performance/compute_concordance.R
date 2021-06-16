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
      Wellcome, wellcome
      
      This R script returns a table with the phenotypes with at least a signfiant result
      that is, p-value < 5e-08 in either the imputed or genotyped data. This script works for 
      Bolt-lmm results and asumes _imputed.stats and _genotype.stats file extention.
      
      Arguments:
      
      --dir= The directory where the summaries are saved
      --ind= The individual 
      --out= The output directory
      
      
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
  stop("Your indivual is missing")
} else {
  ind <- argsL$ind
}
if(is.null(argsL$dir)){
  print.help()
  stop("Your directory is missing")
} else {
  dir <- argsL$dir
}
if(is.null(argsL$out)){
  print.help()
  stop("You should indicate the output directory")
} else {
  out <- argsL$out
}

######################
## DEFINE FUNCTIONS ##
######################
bin_alleles <- function(x){
  if(x <= 1.0 & x > 0.9) return("1.0")
  if(x <= 0.9 & x > 0.8) return("0.9")
  if(x <= 0.8 & x > 0.7) return("0.8")
  if(x <= 0.7 & x > 0.6) return("0.7")
  if(x <= 0.6 & x > 0.5) return("0.6")
  if(x <= 0.5 & x > 0.4) return("0.5")
  if(x <= 0.4 & x > 0.3) return("0.4")
  if(x <= 0.3 & x > 0.2) return('0.3')
  if(x <= 0.2 & x > 0.1) return('0.2')
  if(x <= 0.1 & x > 0.05) return('0.1')
  if(x <= 0.05 & x > 0.04) return('0.05')
  if(x <= 0.04 & x > 0.03) return('0.04')
  if(x <= 0.03 & x > 0.02) return('0.03')
  if(x <= 0.02 & x > 0.01) return('0.02')
  if(x <= 0.01 & x > 0.009) return('0.01')
  if(x <= 0.009 & x > 0.008) return('0.009')
  if(x <= 0.008 & x > 0.007) return('0.008')
  if(x <= 0.007 & x > 0.006) return('0.007')
  if(x <= 0.006 & x > 0.005) return('0.006')
  if(x <= 0.005 & x > 0.003) return('0.005')
  if(x <= 0.003) return('0.003')
}

lvls <- c("1.0","0.9","0.8","0.7","0.6","0.5", "0.4", "0.3", "0.2", "0.1", "0.05", 
          "0.04", "0.03", "0.02", "0.01", "0.009", "0.008", "0.007", "0.006","0.005", "0.003")
ind <- as.numeric(ind)

####################
### START SCRIPT ###
####################

#Read popinfo to assign a population to the individual
amr_ids <- fread("/well/hill/adriank/MEX_BImp/data/array_positions/AMR_ind.txt")
pop <- amr_ids[ind]$POP

#Read the dosage file
in_file <- paste0(dir,"/ind", ind, "_completedata.txt")
dosage_dat <- fread(in_file)

#We are going to make two summaries. The first one will be using the AMR as raw variant frequencies
#[0-1] and the other will be computing MAF [0-0.5]. The difference between this summary and the previous
#is that first we are going to only analyze those variants that had an info score >= 0.3.
#Also, this summary will include rsquared and concordance. 

dosage_dat <- dosage_dat[info >= 0.3]
dosage_dat[, MAF := ifelse(AMR <= 0.5, AMR, 1 - AMR)]
dosage_dat[, `:=` (bin = factor(sapply(MAF, bin_alleles), levels = lvls))]
dosage_dat[, `:=`(bin_amr = factor(sapply(AMR, bin_alleles), levels = lvls))]
dosage_dat[, hard_geno := NA]
dosage_dat[, hard_geno := ifelse(AA >= .9, "AA", hard_geno)]
dosage_dat[, hard_geno := ifelse(AB >= .9, "AB", hard_geno)]
dosage_dat[, hard_geno := ifelse(BB >= .9, "BB", hard_geno)]
dosage_dat[, obs := NA]
dosage_dat[, obs := ifelse(AA_v == 1, "AA", obs)]
dosage_dat[, obs := ifelse(AB_v == 1, "AB", obs)]
dosage_dat[, obs := ifelse(BB_v == 1, "BB", obs)]
dosage_dat[, concordance := (ifelse(hard_geno == obs, 1, 0))]

sum <- dosage_dat[, .(Tot = .N, Conc = sum(concordance, na.rm = T), rho = cor(impute_dosage_kgp, obs_dosage_kgp)), 
                   by = c("bin", "ancestry")]
sum[,concor_per := round(Conc/Tot, 3)]
sum[, Pop := pop]

sum_amr <- dosage_dat[, .(Tot = .N, Conc = sum(concordance, na.rm = T), rho = cor(impute_dosage_kgp, obs_dosage_kgp)), 
                      by = c("bin_amr", "ancestry")]
sum_amr[,concor_per := round(Conc/Tot, 3)]
sum_amr[, Pop := pop]

out_file <- paste0(out, "/", pop,  "/ind_", ind, "_concormaf.txt")
out_fileamr <- paste0(out, "/", pop,  "/ind_", ind, "_concoramr.txt")
fwrite(sum, file = out_file, sep = "\t")
fwrite(sum_amr, file = out_fileamr, sep = "\t")