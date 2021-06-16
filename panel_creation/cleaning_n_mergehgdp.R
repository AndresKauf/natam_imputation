## 
## Este codigo es para hacer una comparacion entre las posiciones de hgdp con las de nuestro panel
## previamente construido. Aqui definire el conjunto de todas las posicones removidas en inmegen12g
## sgdp y hgdp; la union de las posiciones conservadas entre los tres paneles y finalmente definire
## la resta de ambos conjutos para definir las posiciones finales que tendra el panel final. 


library(data.table)

get_chr_pos <- function(data, chr){
  data <- as.data.table(data)
  return(data[CHROM == as.character(chr)][["POS"]])
}
get_chr_pos_hgdp <- function(data, chr){
  data <- as.data.table(data)
  name <- paste0("chr", as.character(chr))
  return(data[CHROM == name][["POS"]])
}

#---- Definiendo la union de todas las posiciones descartadas ----

# Posiciones desechadas en inmegen12g:
inmegen_rfiles <- list.files("./data/genomic_data/inmegen_12g/kept_n_removed/", 
                             pattern = "removed.sites", full.names = T)
inmegen_rdata <- lapply(inmegen_rfiles, fread)
inmegen_removed <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  chr_rm <- lapply(inmegen_rdata, get_chr_pos, i)
  uni <- Reduce(union, chr_rm)
  
  #Make data table with the union of the positions and append to the final data table
  temp_uni <- data.table(chr = i, pos = uni)
  inmegen_removed <- rbindlist(list(inmegen_removed, temp_uni))
}
print(nrow(inmegen_removed))
fwrite(inmegen_removed, file = "./data/genomic_data/inmegen_12g/kept_n_removed/total_removed.txt",
       sep = "\t")

#Posiciones descartadas en sgdp:
sgdp_rfiles <- list.files("./data/genomic_data/simons_americans/kept_n_removed",
                          pattern = "removed.sites", full.names = T)
sgdp_rdata <- lapply(sgdp_rfiles, fread)
sgdp_removed <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  chr_rm <- lapply(sgdp_rdata, get_chr_pos, i)
  uni <- Reduce(union, chr_rm)
  
  #Make data table with the union of the positions and append to the final data table
  temp_uni <- data.table(chr = i, pos = uni)
  sgdp_removed <- rbindlist(list(sgdp_removed, temp_uni))
}
print(nrow(sgdp_removed))
fwrite(sgdp_removed, file = "./data/genomic_data/simons_americans/kept_n_removed/total_removed.txt",
       sep = "\t")

#Posiciones descartadas en hgdp:
hgdp_rfiles <- list.files("./data/genomic_data/HGDP_americans/kept_n_removed/",
                          pattern = "removed.sites", full.names = T)
hgdp_rdata <- lapply(hgdp_rfiles, fread)
hgdp_removed <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  chr_rm <- lapply(hgdp_rdata, get_chr_pos_hgdp, i)
  uni <- Reduce(union, chr_rm)
  
  #Make data table with the union of the positions and append to the final data table
  temp_uni <- data.table(chr = i, pos = uni)
  hgdp_removed <- rbindlist(list(hgdp_removed, temp_uni))
}
print(nrow(hgdp_removed))
fwrite(sgdp_removed, file = "./data/genomic_data/HGDP_americans/kept_n_removed/total_removed.txt",
       sep = "\t")


#Definiendo la union del total de posiciones removidas

total_removed <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  inmegen_chr <- inmegen_removed[chr == i][["pos"]]
  sgdp_chr <- sgdp_removed[chr == i][["pos"]]
  hgdp_chr <- hgdp_removed[chr == i][["pos"]]
  
  sg_mx <- union(inmegen_chr, sgdp_chr)
  sg_mx_hg <- union(sg_mx, hgdp_chr)
  
  temp <- data.table(chr = i, pos = sg_mx_hg)
  total_removed <- rbindlist(list(total_removed, temp))
}
print(nrow(total_removed))
fwrite(total_removed, file = "./data/genomic_data/hg-sg-mex/total_removed.txt",
       sep = "\t")

#---- Definiendo la union de todas las posiciones conservadas ----

#Posiciones conservadas en sgmx
sgmx_pos <- fread("./data/genomic_data/sgdp-mex/sgdp-mex_freeze_keep_positions.txt")

hgdp_kfiles <- list.files("./data/genomic_data/HGDP_americans/kept_n_removed/",
                          pattern = "kept.sites", full.names = T)
hgdp_kdata <- lapply(hgdp_kfiles, fread)
hgdp_kept <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  chr_kept <- lapply(hgdp_kdata, get_chr_pos, i)
  uni <- Reduce(union, chr_kept)
  
  #Make data table with the union of the positions and append to the final data table
  temp_uni <- data.table(chr = i, pos = uni)
  hgdp_kept <- rbindlist(list(hgdp_kept, temp_uni))
}

print(nrow(hgdp_kept))
fwrite(hgdp_kept, file = "./data/genomic_data/HGDP_americans/kept_n_removed/total_kept.txt",
       sep = "\t")

#Definiendo la union de los sitios conservados totales.
total_kept <- data.table(chr = integer(), pos = integer())
for(i in 1:22){
  sgmx_chr <- sgmx_pos[V1 == i][["V2"]]
  hgdp_chr <- hgdp_kept[chr == i][["pos"]]
  uni <- union(sgmx_chr, hgdp_chr)
  temp <- data.table(chr = i, pos = uni)
  total_kept <- rbindlist(list(total_kept, temp))
}

print(nrow(total_kept))


#----Definiendo la resta de los conservados menos los removidos----

positions_t_keep <- data.table(chr = integer(), pos = integer())
positions_t_remove <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  kept_chr <- total_kept[chr ==i][["pos"]]
  rm_chr <- total_removed[chr == i][["pos"]]
  
  t_remove <- intersect(kept_chr, rm_chr)
  t_keep <- setdiff(kept_chr, t_remove)
  
  temp_keep <- data.table(chr = i, pos = t_keep)
  positions_t_keep <- rbindlist(list(positions_t_keep, temp_keep))
  
  temp_remove <- data.table(chr = i, pos = t_remove)
  positions_t_remove <- rbindlist(list(positions_t_remove, temp_remove))
}

print(nrow(positions_t_keep))
print(nrow(positions_t_remove))

for (i in 1:22){
  chr_pos <- positions_t_keep[chr == i]
  name <- paste0("./data/genomic_data/hg-sg-mex/keep_chr", i, ".txt")
  fwrite(chr_pos, file = name, sep = "\t", col.names = F)
}

total_kept 

kgp_pos <- fread("/data/reference_panels/1KGP_plink/ALL.atDNA.biAllelicSNPnoDI.genotypes.id.bim")
kgp_pos[, c("V2", "V3", "V5", "V6") := NULL]


kgp_unique <- data.table(chr = integer(), pos = integer())
hsgmx_unique <- data.table(chr = integer(), pos = integer())
shared <- data.table(chr = integer(), pos = integer())

for (i in 1:22){
  kgp_chr <- kgp_pos[V1 == i][["V4"]]
  hgsmx_chr <- total_kept[chr == i][["pos"]]
  
  kuniqe <- setdiff(kgp_chr, hgsmx_chr)
  hgunique <- setdiff(hgsmx_chr, kgp_chr)
  share <- intersect(kgp_chr, hgsmx_chr)
  
  kunique_dt <- data.table(chr = i, pos = kuniqe)
  kgp_unique <- rbindlist(list(kgp_unique, kunique_dt))
  
  hgunique_dt <- data.table(chr = i, pos = hgunique)
  hsgmx_unique <- rbindlist(list(hsgmx_unique, hgunique_dt))
  
  share_t <- data.table(chr = i, pos = share)
  shared <- rbindlist(list(shared, share_t))
}










