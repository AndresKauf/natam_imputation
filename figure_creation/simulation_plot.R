library(data.table)
library(ggplot2)
library(ggpubr)

color_blind <- c("#009292","#6DB6FF","#CC6677", "#FFB677","#DBD100", "#AA4499") #final palette
color_blindall <- c("#BBBBBB","#332288","#0077BB", "#88CCEE", "#44AA99", "#117733", "#999933",
                    "#DDCC77", "#EE8866","#CC6677", "#882255", "#AA4499")
kauf_theme <- theme(axis.line = element_line(size = .5, colour = "black"), 
                    panel.grid.major = element_line(colour = "#d3d3d3", size = 0.2), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), panel.border = element_blank(),
                    legend.background = element_blank(), legend.box.background = element_blank(),
                    text = element_text(family = "Helvetica"),
                    axis.title.x = element_text(face = "bold",size = 20, color = "black"),
                    axis.title.y = element_text(face = "bold", size = 20, color = "black"),
                    axis.text = element_text(face = "bold", size = 12, color = "black"),
                    legend.title = element_text(face = "bold", size = 20, color ="black" ),
                    legend.text = element_text(size = 18), 
                    strip.text.x = element_text(face = "bold", size = 16, color = "black")) 

xnames <- c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.04", "0.03", "0.02", "0.009",
            "0.008", "0.006", "0.005", "0.003")
xlabs <- c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.04", "0.02",
           "0.005", "0.003")

refs <- paste0("ref", c(0,100,134,200,400,600,800,1000,1500,2000,3000))

nats_n_eur <- fread("~/Documents/moreno_lab/imputation_proj/sim_imp/nats_n_eur.txt")
diffs_natnat <- fread("~/Documents/moreno_lab/imputation_proj/sim_imp/diffs_natnat.txt")

nats_n_eur[, Reference_panel := factor(Reference_panel, levels = refs)]
comp_refs <- refs[2:length(refs)]
diffs_natnat[, Reference_panel := factor(Reference_panel, levels = comp_refs)]

###Compare European vs all nats

compnat_eur <- ggplot(nats_n_eur, aes(x = index, y = mean_rho)) +
  geom_line(aes(linetype = ancestry, color = Reference_panel, size = ancestry))+
  geom_point(aes(group = ancestry)) +
  scale_size_manual(values = c("NAT_NAT" = .8, "EUR_EUR" = 2), guide = "none") +
  scale_x_continuous(breaks = 1:length(xlabs), labels = as.character(xlabs)) +
  scale_y_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9,1)) +
  scale_color_manual(values = color_blindall, name = "Reference") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Ancestry") +
  coord_cartesian(ylim = c(0.5, 1)) +
  ylab( expression(paste("Mean ", r^{2}))) +
  xlab("Variant Frequency") +
  kauf_theme 




natnat_diffplot <- ggplot(diffs_natnat, aes(x = Reference_panel, y = Diff, group = 1)) +
  geom_line(color = "#AA4499") +
  geom_point(color = "#AA4499") +
  ylab("Accuracy increase") +
  xlab("Reference panel") +
  scale_x_discrete( labels = c("100","134","200", "400", "600", "800", "1000", "1500",
                               "2000", "3000")) +
  scale_y_continuous(breaks = seq(from=0, to=0.3, by = 0.02)) +
  kauf_theme +
  facet_wrap(~Var, ncol = 3)

figure3 <- ggarrange(compnat_eur, natnat_diffplot, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave(filename = "~/Documents/moreno_lab/imputation_proj/presentaciones/paper_figures/figure3.hq.jpeg",
       plot = figure3, width = 17, height = 22, units = "cm", 
       dpi = 300, device = "jpeg")
