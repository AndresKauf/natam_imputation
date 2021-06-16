#!/usr/bin/env Rscript

#########################
## Plot and Write filt ##
#########################

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



####################
### START SCRIPT ###
####################

pops <- c("CLM", "MXL", "PEL", "PUR")

for(pop in pops){
  indir <- paste0(dir, "/", pop)
  #If dir is not empty
  if(length(list.files(indir)) != 0){
    #read infomaf files and bind them
    infomaf_files <- list.files(path = indir, pattern = "_infomaf.txt$", full.names = T)
    infomaf <- lapply(infomaf_files, fread)
    infomaf <- rbindlist(infomaf)
    above_maf <- round(nrow(infomaf[MAF >= 0.01])/nrow(infomaf),2)
    above_info <- round(nrow(infomaf[info >= 0.3])/nrow(infomaf),2)
    #Define outnames for plots
    out_maf <- paste0(out, "/", pop, ".MAFdistrib.png")
    out_info <- paste0(out, "/", pop, ".infodistrib.png")
    main_maf <- paste(pop, "imputed MAF distribution")
    main_info <- paste(pop, "imputed info distrubution")
    
    #Plot MAF
    png(out_maf)
    hist(infomaf$MAF, col = "dodgerblue2",
         xlab = "MAF",
         main = main_maf,
         ylim = c(0,1e6),
         breaks = 100)
    abline(v = 0.01, col = "red", lwd = 3, lty = 2)
    legend("topright", legend = paste("MAF > 0.01 =", above_maf),bty ="n", pch=NA)
    dev.off()
    
    #Plot info
    png(out_info)
    hist(infomaf$info, col = "dodgerblue2",
         xlab = "Info value",
         main = main_info)
    abline(v = 0.3, col="firebrick3", lwd=3, lty=2)
    legend("topright", legend = paste("Rsq > 0.3 =", above_info),bty ="n", pch=NA)
    dev.off()
    
    #Now read filter data and write final file
    filterfiles <- list.files(path = indir, pattern = "_filters.txt$", full.names = T)
    filters <- lapply(filterfiles, fread)
    filters <- rbindlist(filters)
    promedios <- filters[, .(maf_filt = round(mean(MAF)), 
                             info_filt = round(mean(INFO)), all_filt = round(mean(ALL)))]
    fwrite(promedios, file = paste0(indir, "/all_filter_sum.txt"), sep = "\t")
    
  }else{
    print("No data in dir")
  }
}