library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggpubr")
library("data.table")
library("tidyr")
library("cowplot")

##### Plotting map ----
world <- ne_countries(scale = "medium", returnclass = "sf")
sites_sgdp <- data.frame(latitude = c(-29.53, -8.40, -9.79, 5.26, 17.5, 20.46,
                                       30.44, -7.01, 16.96 ), 
                         longitude = c(-63.1, -53.61, -64.33, -73.3, -99.05, -88.94,
                                      -111.26, -76.64, -94.96))
sites_hgdp <- data.frame(latitude = c(18.67, 4.42, -6.76, -8.68, 30.51, 20.31,
                                       28.72), 
                         longitude = c(-90.15, -71.89, -56.07, -66.97, -110.12, -98.30,
                                      -107.51))
sites_inmegen <- data.frame(latitude = c(19.50, 20.45, 25.60, 16.12), 
                            longitude = c(-89.10, -99.66, -107.33, -96.52))
sites_mxb <- data.frame(latitude = c(23.54, 17.87, 22, 19.70, 25.93, 18.28), 
                        longitude = c(-102.93, -92.64, -102.40, -97.39, -99.23,
                                     -101.08))

map <- ggplot(data = world) +
  geom_sf(fill = "#F7F7F7") +
  geom_point(data = sites_hgdp, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#88CCEE") +
  geom_point(data = sites_sgdp, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#117733")+
  geom_point(data = sites_inmegen, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#DDCC77")+
  geom_point(data = sites_mxb, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#CC6677")+
  annotate(geom = "point", size = 5, colour = "#88CCEE", x = -114, y = 5) +
  annotate(geom = "text", label = "HGDP", x = -103, y = 5, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#117733", x = -114, y = 0) +
  annotate(geom = "text", label = "SGDP", x = -103, y = 0, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#DDCC77", x = -114, y = -5) +
  annotate(geom = "text", label = "INMEGEN", x = -97.5, y = -5, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#CC6677", x = -114, y = -10) +
  annotate(geom = "text", label = "Mexico", x = -101.5, y = -10, fontface = "bold", size = 5) +
  coord_sf(xlim = c(-118.53, -33.45), ylim = c(-35.35, 32.50), expand = F) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "#FFFFFF"), axis.text = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank()) 

map2 <- ggplot(data = world) +
  geom_sf(fill = "#F7F7F7") +
  geom_point(data = sites_hgdp, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#88CCEE") +
  geom_point(data = sites_sgdp, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#117733")+
  geom_point(data = sites_inmegen, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#DDCC77")+
  geom_point(data = sites_mxb, aes(x = longitude, y = latitude),
             size = 3, shape = 19,  col = "#CC6677")+
  annotate(geom = "point", size = 5, colour = "#88CCEE", x = -114, y = -10) +
  annotate(geom = "text", label = "HGDP", x = -108, y = -10, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#117733", x = -114, y = -15) +
  annotate(geom = "text", label = "SGDP", x = -108, y = -15, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#DDCC77", x = -114, y = -20) +
  annotate(geom = "text", label = "INMEGEN", x = -106, y = -20, fontface = "bold", size = 5) +
  annotate(geom = "point", size = 5, colour = "#CC6677", x = -114, y = -25) +
  annotate(geom = "text", label = "MX Biobank", x = -105, y = -25, fontface = "bold", size = 5) +
  coord_sf(xlim = c(-118.53, -33.45), ylim = c(-35.35, 32.50), expand = F) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "#FFFFFF"), axis.text = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank()) 



##### Plotting props ----

data <- fread("~/Documents/moreno_lab/imputation_proj/data/data_forfigs/summary_table.txt")

data <- data[23,]
data_plt <- data[, .(kgp_unq = kgp_unique/total, inter = intersection/total, natpan_unq = natpan_unique/total)]
data_plt[, index := 1]
data_plt <- gather(data_plt, type, proportion, -index)
data_plt <- as.data.table(data_plt)
data_plt <- data_plt[order(proportion, decreasing = T)]
data_plt[, type := factor(type, levels = c("natpan_unq", "inter", "kgp_unq"))]

colorsprop <- c("#AA4499", "#DDDDDD", "#555555")
colorsprop2 <- c("#DDAA33", "#BB5566", "#004488")
prop_plot <- ggplot(data_plt, aes(x=index, y = proportion, fill = type)) +
  geom_bar(stat = "identity", width = .5) +
  scale_fill_manual(values = colorsprop2, name = "Type of SNP", 
                    labels = c("NATS unique", "Intersection", "1KGP unique") ) +
  scale_x_continuous(breaks = c(1), labels = c("1KGP + NATS")) +
  scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),expand = c(0,0), limits = c(0,1.01)) +
  theme_bw() +
  labs(x="", y = "Proportion") +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
        legend.position = "right", legend.direction = "vertical",
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(face = "bold", size = 18, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 18, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 18)) 

##### Plotting admix ----
#Read new file
q_data <- fread("~/Documents/moreno_lab/imputation_proj/data/data_forfigs/natpan.cyamr.clean.tidydat.3.Q")
#Set index
q_data[, index:=seq(1,nrow(q_data))]
q_data[, id := NULL]
#gather
q_data_new <- gather(q_data, admix_group, porportion, -pop, -index)


#Plotting

colors <- c("#009292", "#FFB677", "#AA4499")
admix_plot <- ggplot(q_data_new, aes(x=index, y = porportion, fill = admix_group )) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors, name = "Genetic Ancestry") +
  scale_x_continuous(breaks = c(73,200,308,540), labels = c("NAT", "EUR", "AFR", "AMR")) +
  scale_y_continuous(breaks = NULL,expand = c(0,0), limits = c(0,1.01)) +
  theme_bw() +
  labs(x="", y = "") +
  coord_cartesian(xlim=c(-1,700), ylim = c(-0.001, 1)) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"), axis.line.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(face = "bold", size = 6, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 10)) +
  annotate(x=73, xend = 540, y = 0, yend = 0, colour = "black", lwd = 0.5, geom = "segment") +
  theme(legend.position = "none") 

admix_plot_comp <- ggplot(q_data_new, aes(x=index, y = porportion, fill = admix_group )) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors, name = "Genetic Ancestry", labels = c("African", "European",
                                                                           "Native American")) +
  scale_x_continuous(breaks = c(73,200,308,540), labels = c("NATS", "EUR", "AFR", "AMR")) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1),expand = c(0,0), limits = c(0,1.01)) +
  theme_bw() +
  labs(x="Population", y = "Ancestry Proportion") +
  coord_cartesian(xlim=c(-1,700), ylim = c(-0.001, 1)) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_line(colour = "black"), axis.line.y = element_line(colour="black"),
        legend.position = "bottom", legend.direction = "horizontal",
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(face = "bold", size = 18, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 18, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 20)) +
  annotate(x=73, xend = 540, y = 0, yend = 0, colour = "black", lwd = 0.5, geom = "segment") +
  theme(legend.position = "bottom") 

#####Plotting them together ----

figure1 <- ggdraw() +
  draw_plot(map2, x = 0, y = 0.5, width = .65, height = .5) +
  draw_plot(prop_plot, x = 0.65, y = 0.5, width = .35, height = .5) +
  draw_plot(admix_plot, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

ggsave(filename = "~/Documents/moreno_lab/imputation_proj/presentaciones/paper_figures/figure1.hq.jpeg",
       plot = figure1, width = 16, height = 15, units = "cm", 
       dpi = 300, device = "jpeg")
