library(data.table)
library(ggplot2)
library(ggpubr)

amr_info <- fread("~/Documents/moreno_lab/imputation_proj/data/amr_infothresh.txt")
#Statistical tests for each population
mxl <- amr_info[POP == "MXL"]
mtest <- t.test(mxl$thresh_one, mxl$thresh_two, paired = T)$p.value
pel <- amr_info[POP == "PEL"]
ptest <- t.test(pel$thresh_one, pel$thresh_two, paired = T)$p.value
clm <- amr_info[POP == "CLM"]
clmtest <- t.test(clm$thresh_one, clm$thresh_two, paired = T)$p.value
pur <- amr_info[POP == "PUR"]
purtest <- t.test(pur$thresh_one, pur$thresh_two, paired = T)$p.value

sp <- ggscatter(amr_info, x = "NAT", y = "diff",
                add = "reg.line",
                conf.int = TRUE,                # Add confidence interval
                palette = "jco",
                color = "#225555",
                ylab = "Increase of SNPs above quality threshold",
                xlab = "Proportion of Native American ancestry") +
    stat_cor( label.x = 0.1)
  xlab("Increase of SNPs above quality threshold") +
  ylab("Proportion of Native American ancestry")

  
real_sfs <- fread("~/Documents/moreno_lab/imputation_proj/data/real_sfs.txt")
sim_sfs <- fread("~/Documents/moreno_lab/imputation_proj/data/sim_sfs.txt")
sim_mean <- 2173.766
mega_meand <- 2308.403
binf <- c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.04", "0.03", "0.02", "0.01",
          "0.009", "0.008", "0.007", "0.006", "0.005", "0.003")
sim_sfs[, bin := factor(bin, levels=binf)]
sim_sfs <- sim_sfs[order(bin)]

real_sfs[, bin := factor(bin, levels = binf)]
real_sfs <- real_sfs[order(real_sfs)]

#Plot real and simulated sfs
kauf_theme <- theme(axis.line = element_line(size = .5, colour = "black"), 
                    panel.grid.major = element_line(colour = "#d3d3d3", size = 0.2), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), panel.border = element_blank(),
                    legend.background = element_blank(), legend.box.background = element_blank(),
                    text = element_text(family = "Helvetica"),
                    axis.title.x = element_text(face = "bold",size = 16, color = "black"),
                    axis.title.y = element_text(face = "bold", size = 16, color = "black"),
                    axis.text = element_text(face = "bold", size = 12, color = "black"),
                    legend.title = element_text(face = "bold", size = 24, color ="black" ),
                    legend.text = element_text(face = "bold", size = 18), 
                    strip.text.x = element_text(face = "bold", size = 16, color = "black")) 

mega_meand <- round(mega_meand, 0)
plot_real <- ggplot(real_sfs, aes(x = bin, y = prop)) +
  geom_bar(stat = "identity", fill = "#009292") +
  labs(x = "Minor Allele Frequency", y = "Proportion") +
  annotate("text", x = "0.1", y = 0.3, label = paste("Mean distance between markers =", mega_meand)) +
  ggtitle("MEGA Frequency Spectrum for chromosome 9") +
  kauf_theme

sim_mean <- round(sim_mean,0)
plot_sim <- ggplot(sim_sfs, aes(x = bin, y = prop)) +
  geom_bar(stat = "identity", fill = "#FFB677") +
  labs(x = "Minor Allele Frequency", y = "Proportion") +
  annotate("text", x = "0.1", y = 0.3, label = paste("Mean distance between markers =", sim_mean)) +
  ggtitle("Simulated Frequency Spectrum for chromosome 9") +
  kauf_theme

ggarrange(plot_real, plot_sim,
          ncol = 1, nrow = 2,
          labels = c("A", "B"))
