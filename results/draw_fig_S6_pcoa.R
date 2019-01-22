###############################
#
# Figure S5
# PCoA of communities
#
###############################

library(tidyr)
library(ggplot2)
library(cowplot)


# read in data
pcoa.ca <- read.table(file="data/mothur/phage_ecol.0.03.ceca_exp.braycurtis.0.03.lt.pcoa.axes", 
                      sep="\t", header=TRUE)[1:3]

pcoa.il <- read.table(file="data/mothur/phage_ecol.0.03.il_exp.braycurtis.0.03.lt.pcoa.axes", 
                      sep="\t", header=TRUE)[1:3]

# format for metadata
pcoa.ca.meta <- separate(pcoa.ca, group, 
                         c("site", "treatment", "dpi", "replicate"), "_")
pcoa.il.meta <- separate(pcoa.il, group, 
                         c("site", "treatment", "dpi", "replicate"), "_")


p.ceca.pcoa <- ggplot(pcoa.ca.meta, aes(x=axis1, y=axis2, 
                      color=factor(treatment, labels=c("Cj", "Cj_phg", "Control")), 
                       shape=factor(dpi, labels=c("1", "2", "3", "4", "5")))) +
                  geom_point() +
                  labs(color="group", shape="dpt") +
                  ylab("PCoA 2 (22.4%)") +
                  xlab("PCoA 1 (10%)") +
                  scale_shape_manual(values=c(8, 15, 16, 17, 18)) +
                  theme_bw() +
                  theme(panel.grid.major=element_blank(), 
                        panel.grid.minor=element_blank(),
                        legend.position=c(0.345, 0.914), 
                        legend.background=element_blank(),
                        legend.text=element_text(size=8), legend.title=element_text(size=8),
                        legend.box.background=element_rect(colour="black", size = 0.3),
                        legend.direction="horizontal",
                        legend.spacing.y=unit(0.01, "mm"),
                        legend.margin=margin(c(0.5,2,0.5,1), unit="mm")) +
                  scale_x_continuous(breaks=seq(-0.5, 0.5, 0.25), limits=c(-0.5,0.5), 
                     expand=c(0,0)) +
                  scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3), limits=c(-0.6, 0.6), 
                     expand=c(0,0))

p.ileum.pcoa <- ggplot(pcoa.il.meta, aes(x=axis1, y=axis2, 
                        color=factor(treatment, labels=c("Cj", "Cj_phg", "Control")), 
                        shape=factor(dpi, labels=c("2", "3", "4", "5")))) +
                    geom_point() +
                    labs(color="group", shape="dpt") +
                    ylab("PCoA 2 (13.1%)") +
                    xlab("PCoA 1 (15.9%)") +
                    scale_shape_manual(values=c(15, 16, 17, 18, 19)) +
                    theme_bw() +
                    theme(panel.grid.major=element_blank(), 
                      panel.grid.minor=element_blank(), 
                      legend.position=c(0.336, 0.914), 
                       legend.background=element_blank(),
                        legend.text=element_text(size=8), legend.title=element_text(size=8),
                         legend.box.background=element_rect(colour = "black", size = 0.3),
                          legend.direction="horizontal",
                           legend.spacing.y = unit(0.01, "mm"),
                            legend.margin=margin(c(0.5,2,0.5,1), unit="mm")) +
                    scale_x_continuous(breaks=seq(-0.5, 0.5, 0.25), limits=c(-0.5,0.5), 
                     expand=c(0,0)) +
                    scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3), limits=c(-0.6, 0.6), 
                     expand=c(0,0))




tiff(filename="results/figures/figS6_pcoa_fix.tiff", 
     width=210, height=100, units="mm", res=300)
plot_grid(p.ileum.pcoa, p.ceca.pcoa, labels = c('A', 'B'), scale = 0.97)
dev.off()