###############################
#
# Figure S5
# PCoA of communities
#
###############################

library("tidyr")
library("ggplot2")


# read in data
pcoa.ca <- read.table(file="data/mothur/phage_ecol.0.03.ceca_exp.braycurtis.0.03.lt.pcoa.axes", sep="\t", header=TRUE)[1:3]

pcoa.il <- read.table(file="data/mothur/phage_ecol.0.03.il_exp.braycurtis.0.03.lt.pcoa.axes", sep="\t", header=TRUE)[1:3]

# format for metadata
pcoa.ca.meta <- separate(pcoa.ca, group, c("site", "treatment", "dpi", "replicate"), "_")
pcoa.il.meta <- separate(pcoa.il, group, c("site", "treatment", "dpi", "replicate"), "_")


p.ceca.pcoa <- ggplot(pcoa.ca.meta, aes(x=axis1, y=axis2, 
                        color = factor(treatment, labels=c("Cj", "Cj_phg", "Control")), 
                        shape = factor(dpi, labels=c("1", "2", "3", "4", "5")))) +
                geom_point() +
                labs(color = "group", shape = "dpt") +
                ylab("PCoA 2 (22.4%)") +
                xlab("PCoA 1 (10%)") +
                scale_shape_manual(values=c(15, 16, 17, 18, 19)) +
                theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p.ileum.pcoa <- ggplot(pcoa.il.meta, aes(x=axis1, y=axis2, 
                        color = factor(treatment, labels=c("Cj", "Cj_phg", "Control")), 
                        shape = factor(dpi, labels=c("2", "3", "4", "5")))) +
                  geom_point() +
                  labs(color = "group", shape = "dpt") +
                  ylab("PCoA 2 (13.1%)") +
                  xlab("PCoA 1 (15.9%)") +
                  scale_shape_manual(values=c(15, 16, 17, 18, 19)) +
                  theme_bw() +
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




tiff(filename="results/figures/figS6a_ileum_pcoa.tiff", width=140, height=100, units="mm", res=300)
p.ileum.pcoa
dev.off()

tiff(filename="results/figures/figS6b_ceca_pcoa.tiff", width=140, height=100, units="mm", res=300)
p.ceca.pcoa
dev.off()