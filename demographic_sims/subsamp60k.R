library(data.table)
"%nin%" <- Negate("%in%")
pos <- fread("~/projects/sim_imp2000/pos.maf0.01.txt")
pos <- pos$V1
dups <- pos[duplicated(pos)]
pos <- pos[pos %nin% dups]
samp_pos <- sample(pos, 60000)
samp_pos <- sort(samp_pos)
new <- data.table(chr = 9, pos = samp_pos)
fwrite(new, file = "~/projects/sim_imp2000/local_ancestry/data/sampled60kpos.txt", sep = "\t",
       col.names = F)

admix <- fread("~/projects/sim_imp/local_ancestry2/data/admix.target.txt")
a1 <- admix[1:100,]
a2 <- admix[101:200,]
a3 <- admix[201:300,]

fwrite(a1, file = "~/projects/sim_imp2000/local_ancestry/data/admix.target1.txt",
       col.names = F, sep = "\t")
fwrite(a2, file = "~/projects/sim_imp2000/local_ancestry/data/admix.target2.txt",
       col.names = F, sep = "\t")
fwrite(a3, file = "~/projects/sim_imp2000/local_ancestry/data/admix.target3.txt",
       col.names = F, sep = "\t")

col <- fread("~/projects/sim_imp/local_ancestry2/data/admix.col.txt", header = F)
a1 <- col[1:100,]
a2 <- col[101:200,]
a3 <- col[201:300,]

fwrite(a1, file = "~/projects/sim_imp2000/local_ancestry/data/admix.col1.txt",
       col.names = F, sep = "\t")
fwrite(a2, file = "~/projects/sim_imp2000/local_ancestry/data/admix.col2.txt",
       col.names = F, sep = "\t")
fwrite(a3, file = "~/projects/sim_imp2000/local_ancestry/data/admix.col3.txt",
       col.names = F, sep = "\t")
