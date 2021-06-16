library(data.table)

frqx <- fread("~/projects/sim_imp2000/imputation_run/freq_files/admix.frqx")
frqx <- frqx[,1:7]
setnames(frqx, colnames(frqx), c("chr", "position", "A1", "A2", "homA1", "het", "homA2"))
frqx[, Var_freq := (2*homA1 + het)/694]
af <- frqx[, c("position", "Var_freq")]
fwrite(af, file = "~/projects/sim_imp2000/imputation_run/freq_files/af.txt", sep = "\t")
