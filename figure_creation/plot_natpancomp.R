library(data.table)
library(ggplot2)
library(tidyr)
library(ggpubr)

"%nin%" <- Negate("%in%")

kgp_dir <- "~/Documents/moreno_lab/imputation_proj/data/hgmx_second/sums_one/"

color2 <- c("tan2", "darkcyan", "dodgerblue2")
colors <- c("cyan", "teal", "olive", "green", "rose", "wine")
color2_hex <-c("#999933", "#117733","#CC6677", "#882255")
colors_hex <- c("#88CCEE", "#44AA99", "#999933", "#117733", "#CC6677", "#882255")
color_blind <- c("#009292","#6DB6FF","#CC6677", "#FFB677","#DBD100", "#AA4499") #final palette
color_blind2 <- c("#FFB677","#AA4499","#DBD100", "#CC6677") #pastel
color_blind3 <- c("#AA4499", "#CC6677")

kauf_theme <- theme(axis.line = element_line(size = .5, colour = "black"), 
                    panel.grid.major = element_line(colour = "#d3d3d3", size = 0.2), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), panel.border = element_blank(),
                    legend.background = element_blank(), legend.box.background = element_blank(),
                    text = element_text(family = "Helvetica"),
                    axis.title.x = element_text(face = "bold",size = 10, color = "black"),
                    axis.title.y = element_text(face = "bold", size = 10, color = "black"),
                    axis.text = element_text(face = "bold", size = 6, color = "black"),
                    legend.title = element_text(face = "bold", size = 8, color ="black" ),
                    legend.text = element_text(size = 8), 
                    strip.text.x = element_text(face = "bold", size = 6, color = "black")) 

xnames <- c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.04", "0.03", "0.02", "0.009",
            "0.008", "0.006", "0.005", "0.003")

kgp_chr2 <- paste0(kgp_dir, "chr2")
kgp_chr9 <- paste0(kgp_dir, "chr9")
kgp_rhocom <- list.files(path = kgp_chr9, pattern = "common_ann.txt", full.names = T)
#kgp_rhocom <- c(list.files(path=kgp_chr2, pattern = "common_ann.txt", full.names = T),
 #              list.files(path=kgp_chr9, pattern = "common_ann.txt", full.names = T))

kgp_rhocom <- lapply(kgp_rhocom, fread)
kgp_rhocom <- rbindlist(kgp_rhocom)
kgp_rhocommeans <- kgp_rhocom[, .(mean_rho = mean(rho, na.rm =T)),
                              by = c("bin", "ancestry", "Pop")]
kgp_rhocommeans <- kgp_rhocommeans[order(ancestry, -bin)]
#Change population names
kgp_rhocommeans[, Pop2 := "NA"]
kgp_rhocommeans[, Pop2 := ifelse(Pop == "MXL", "Mexico", Pop2)]
kgp_rhocommeans[, Pop2 := ifelse(Pop == "CLM", "Colombia", Pop2)]
kgp_rhocommeans[, Pop2 := ifelse(Pop == "PEL", "Peru", Pop2)]
kgp_rhocommeans[, Pop2 := ifelse(Pop == "PUR", "Puerto Rico", Pop2)]
kgp_rhocommeans[, Pop2 := factor(Pop2, levels = c("Peru", "Mexico", "Colombia", "Puerto Rico"))]
kgp_rhocommeans[, Reference_panel := "KGP"]

#Ploting
xlabs <- c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.04", "0.02",
           "0.005", "0.003")
kgp_rhocommeans_c <- kgp_rhocommeans[bin %in% xlabs]
index <- rep(1:10, 6, each = 4)
kgp_rhocommeans_c[,index := index]


perform_kgprhocom <- ggplot(kgp_rhocommeans_c, aes(x = index, y = mean_rho, color = ancestry, group =ancestry)) +
  geom_line(size = 1) +
  geom_point(size = 1.3) +
  scale_x_continuous(breaks = 1:length(xlabs), labels = xlabs) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  ylab( expression(paste("Mean ", r^{2}))) + 
  xlab("Variant Frequency") +
  scale_color_manual(values = color_blind, name = "Ancestry") +
  kauf_theme +
  facet_wrap(~Pop2, ncol = 2)


#Comparing refs
nat_dir <- "~/Documents/moreno_lab/imputation_proj/data/hgmx_second/sums_natpan/"
nat_chr2 <- paste0(nat_dir, "chr2")
nat_chr9 <- paste0(nat_dir, "chr9")
nat_rhocom <- list.files(path = nat_chr9, pattern = "common_ann.txt", full.names = T)
#nat_rhocom <- c(list.files(path=nat_chr2, pattern = "common_ann.txt", full.names = T),
 #               list.files(path=nat_chr9, pattern = "common_ann.txt", full.names = T))
nat_rhocom <- lapply(nat_rhocom, fread)
nat_rhocom <- rbindlist(nat_rhocom)
nat_rhocommeans <- nat_rhocom[, .(mean_rho = mean(rho, na.rm =T)),
                              by = c("bin", "ancestry", "Pop")]
nat_rhocommeans <- nat_rhocommeans[order(ancestry, -bin)]
#Change population names
nat_rhocommeans[, Pop2 := "NA"]
nat_rhocommeans[, Pop2 := ifelse(Pop == "MXL", "Mexico", Pop2)]
nat_rhocommeans[, Pop2 := ifelse(Pop == "CLM", "Colombia", Pop2)]
nat_rhocommeans[, Pop2 := ifelse(Pop == "PEL", "Peru", Pop2)]
nat_rhocommeans[, Pop2 := ifelse(Pop == "PUR", "Puerto Rico", Pop2)]
nat_rhocommeans[, Pop2 := factor(Pop2, levels = c("Peru", "Mexico", "Colombia", "Puerto Rico"))]
nat_rhocommeans[, Reference_panel := "KGP + NATS"]


xlabs2 <- xnames[9:14]
kgp_rhocommeansrare <- kgp_rhocommeans[bin %in% xlabs2]
nat_rhocommeansrare <- nat_rhocommeans[bin %in% xlabs2]
index2 <- rep(1:6,6, each = 4)
kgp_rhocommeansrare[,index := index2]
nat_rhocommeansrare[, index := index2]

kgp_rhocommeansrare <- kgp_rhocommeansrare[ancestry == "NAT_NAT" | ancestry == "EUR_EUR"]
nat_rhocommeansrare <- nat_rhocommeansrare[ancestry == "NAT_NAT" | ancestry == "EUR_EUR"]

kgp_natrhocomrare <- rbindlist(list(kgp_rhocommeansrare, nat_rhocommeansrare))
#kgp_natrhocomrare <- kgp_natrhocomrare[ancestry == "NAT_NAT" | ancestry ==  "EUR_NAT" | ancestry == "AFR_NAT" | ancestry == "EUR_EUR"]
compare_kgprhorare <- ggplot(kgp_natrhocomrare, aes(x = index, y = mean_rho)) +
  geom_line(aes(linetype = Reference_panel, color = ancestry), size = .8)+
  geom_point(aes(group = ancestry)) +
  scale_x_continuous(breaks = 1:length(xlabs2), labels = as.character(xlabs2)) +
  scale_y_continuous(breaks = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  scale_color_manual(values = color_blind2, name = "Ancestry") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Reference Panel") +
  coord_cartesian(ylim = c(0.3, 1)) +
  ylab( expression(paste("Mean ", r^{2}))) +
  xlab("Variant Frequency") +
  kauf_theme +
  theme(legend.key.size =  unit(.8, "cm")) +
  facet_wrap(~Pop2, ncol = 2)

#png(filename = "~/Documents/moreno_lab/imputation_proj/presentaciones/paper_figures/figure2.hq.png",
   # width = 8, height = 12, units = "cm", res = 300)
figure2 <- ggarrange(perform_kgprhocom, compare_kgprhorare, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
ggsave(filename = "~/Documents/moreno_lab/imputation_proj/presentaciones/paper_figures/figure2.hq.jpeg",
       plot = figure2, width = 16, height = 15, units = "cm", 
       dpi = 300, device = "jpeg")

##### Supplementary figure with all ancestries ----

nat_rhocommeans_c <- nat_rhocommeans[bin %in% xlabs]
nat_rhocommeans_c[, index := index]
kgp_natrhocomm <- rbindlist(list(kgp_rhocommeans_c, nat_rhocommeans_c))

compare_kgprhocom <- ggplot(kgp_natrhocomm, aes(x = index, y = mean_rho)) +
  geom_line(aes(linetype = Reference_panel, color = ancestry), size = .8)+
  geom_point(aes(group = ancestry)) +
  scale_x_continuous(breaks = 1:length(xlabs), labels = as.character(xlabs)) +
  scale_y_continuous(breaks = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  scale_color_manual(values = color_blind, name = "Ancestry") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Reference Panel") +
  coord_cartesian(ylim = c(0.3, 1)) +
  ylab( expression(paste("Mean ", r^{2}))) +
  xlab("Variant Frequency") +
  kauf_theme +
  theme(legend.key.size =  unit(2, "cm")) +
  facet_wrap(~Pop2, ncol = 2) 



###### Statistical tests ----

#paste rho from kgp_rhocom to nat_rhocom (they are in the same order)

nat_rhocom[, rho_kgp := kgp_rhocom$rho]
signif <- nat_rhocom[, .(test = t.test(rho, rho_kgp, paired = T)$p.value),
                     by = c("bin", "ancestry","Pop")]
signif <- signif[bin %in% xlabs2]

