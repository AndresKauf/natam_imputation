#!/usr/bin/env Rscript

######################
## Annotate Results ##
######################

#Author: Andres Kaufmann


library(data.table)
library(parallel)

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
      --array = The array of the simulation
      --chr = The chromosome of the simulation
      --val = The dir with the validation data
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
  stop("Your directory is missing")
} else {
  ind <- argsL$ind
}
if(is.null(argsL$dir)){
  print.help()
  stop("You have to provide a directory with the results")
} else{
  dir <- argsL$dir
}
if(is.null(argsL$chr)){
  print.help()
  stop("You have to provide the chromosome")
} else{
  chr <- argsL$chr
}
if(is.null(argsL$array)){
  print.help()
  stop("You have to provide the simulated array")
} else{
  array <- argsL$array
}
if(is.null(argsL$val)){
  print.help()
  stop("You have to provide the validation dir")
} else{
  val <- argsL$val
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
"%nin%" <- Negate("%in%")
#Defining directories and variables
ancestry_dir <- "/well/hill/adriank/MEX_BImp/data/hap_ancestry/haps/"
freq_dir <- "/well/hill/adriank/MEX_BImp/data/1KGP/allele_freqs/"
val_dir <- val
array_dir <- "/well/hill/adriank/MEX_BImp/data/array_positions"
array_file <- paste0(array_dir, "/", array, "_", chr, "_positions.txt")
impute_dir <- paste0(dir, "/", "indiv", ind)
amr_ids <- fread("/well/hill/adriank/MEX_BImp/data/array_positions/AMR_ind.txt")


#Define the AMR id and read the .bed files
id <- amr_ids[ind]$ID
pop <- amr_ids[ind]$POP
a_file <- fread(paste0(ancestry_dir, id, "_A.bed"))
b_file <- fread(paste0(ancestry_dir, id, "_B.bed"))
#Get the data from the chromosome 

a_file <- a_file[V1 == chr]
b_file <- b_file[V1 == chr]

#Cargar las frecuencias alelicas de cada poblacion calculadas previamente

#allele_file <- paste0(freq_dir, pop, "_MAF_", chr, ".txt" )
allele_file <- paste0(freq_dir, "AMR_MAF_", chr, ".txt" )
freqs <- fread(allele_file)
setkey(freqs, "position")
pos <- freqs$position

#Leer el archivo de validacion 
val_file <- paste0(val_dir, "/" , "ind_", ind, "_", array,"_chr", chr,"_validation.gz")
validation <- fread(paste("zcat", val_file))
validation <- validation[,-c("V1", "V2")]
setnames(validation, names(validation), c("position","A_v", "B_v" ,"AA_v", "AB_v", "BB_v"))

#Necesito leer el archivo de dosages y el de info para hacer uno solo que contenga dosages e info
#Despues tengo que unirlo con las frecuencias alelicas AMR de 1kgp

dosage_files <- list.files(path = impute_dir, pattern = "chunk\\d+$", full.names = T)

dosage_dat <- lapply(dosage_files, fread)
dosage_dat <- rbindlist(dosage_dat)
dosage_dat <- dosage_dat[order(V3)]
dosage_dat <- dosage_dat[V3 %in% pos] #posiciones de cada poblacion
dosage_dat <- dosage_dat[,-c("V1")] #Eliminar columna V1
setnames(dosage_dat, names(dosage_dat), c("id", "position", "A", "B", "AA", "AB", "BB"))
#Merge con 1kgp
setkey(dosage_dat, "position")
dosage_dat <- dosage_dat[freqs, nomatch = 0]
dosage_dat[, i.id := NULL]

info_files <- list.files(path = impute_dir, pattern = "_info$", full.names = T)
info_dat <- lapply(info_files, fread)
info_dat <- rbindlist(info_dat)
info_dat <- info_dat[, c("position", "info")]
setkey(info_dat, "position")

#Merge
dosage_dat <- dosage_dat[info_dat, nomatch = 0]

#Leemos las posciones de MEGA y las quitamos de dosage_dat
array_pos <- fread(array_file)
array_pos <- array_pos$V2
dosage_dat <- dosage_dat[position %nin% array_pos]

#Hacer el merge con los datos de validacion
setkey(validation, "position")
setkey(dosage_dat, "position")
#Otro merge
dosage_dat <- dosage_dat[validation, nomatch = 0]
#quitar posiciones que no concuerden
dosage_dat <- dosage_dat[, same := ifelse(A == A_v & B == B_v, TRUE, FALSE)]
dosage_dat <- dosage_dat[same == TRUE] #Datos limpios 
#Drop some columns
dosage_dat[, A_v := NULL]
dosage_dat[, B_v := NULL]
dosage_dat[, same := NULL]
dosage_dat[, a0 := NULL]
dosage_dat[, a1 := NULL]
#Generar funcion que determine la ancestria en A y B de una posicion
#Paralelizar esta cagada:

no_cores <- 10
cl <- makeCluster(no_cores)

ancestry <- function(pos,bed){
  n <- nrow(bed)
  for (i in seq(1,n)){
    if( pos >= bed[i,2] & pos <= bed[i,3] ){
      return(bed[i,4])
    }
  }
  return("UNK")
}

positions <- dosage_dat$position
hap_a <- parSapply(cl,positions, ancestry, a_file)
hap_b <- parSapply(cl, positions, ancestry, b_file)

stopCluster(cl)

dosage_dat[, ancestry_a := hap_a]
dosage_dat[, ancestry_b := hap_b]

dosage_dat[, impute_dosage_kgp := (AB + 2*(BB))]
dosage_dat[, obs_dosage_kgp := (AB_v + 2*(BB_v))]
dosage_dat <- dosage_dat[ancestry_a != "UNK" & ancestry_b != "UNK"]


dosage_dat <- dosage_dat[, ancestry := paste0(ancestry_a, "_",ancestry_b)]
dosage_dat <- dosage_dat[, ancestry := (ifelse(ancestry == "EUR_AFR",
                                               "AFR_EUR",
                                               ancestry))]
dosage_dat <- dosage_dat[, ancestry := (ifelse(ancestry == "NAT_AFR",
                                               "AFR_NAT",
                                               ancestry))]
dosage_dat <- dosage_dat[, ancestry := (ifelse(ancestry == "NAT_EUR",
                                               "EUR_NAT",
                                               ancestry))]

########################
####Write dosage data###
########################

fwrite(dosage_dat, file = paste0(out,"/ind",ind, "_completedata.txt"), sep = "\t")
