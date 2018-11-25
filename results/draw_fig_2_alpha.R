###############################
#
# Figure 2
# Alpha diversity.
#
###############################

library("readxl")
library("ggpubr")
library("cowplot")

## read in data
alpha.div <- read.table(file="data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.groups.summary", sep="\t", header=TRUE) 

meta <- read_excel("data/references/phage_ecol.MIMARKS.survey.host-associated.4.0.xlsx", range="A12:U108")

# clean data
meta.clean <- meta[meta$host_disease!="none",]
meta.clean <- meta.clean[-c(4,29),] # communities "il_cmp_d2_2" and "il_cmpphg_d4_3" removed 
                                    # because of insufficiaent sampling effort

# merge dataframes
alpha.meta <- merge(alpha.div, meta.clean, by.x="group", by.y="*sample_name")
alpha.meta$host_age <- as.numeric(alpha.meta$host_age)

# add dpi
alpha.meta$dpi <- as.factor(alpha.meta$host_age-20)

# add group
alpha.meta$group2 <- as.factor(ifelse(alpha.meta$chem_administration=="mock", "Cj", "Cj_phg"))
# subset
alpha.meta.caeca <- alpha.meta[alpha.meta$`*env_feature`=="CAECUM",]
alpha.meta.ileum <- alpha.meta[alpha.meta$`*env_feature`=="ILEUM",]

## caecum ##

# diversity
# Inverse simpson plots
ca.invsim.p <- ggboxplot(alpha.meta.caeca, x="dpi", y="invsimpson", 
                  fill="group2", palette=c("deepskyblue3", "white")) + 
                theme(legend.position="none", axis.title.x=element_blank()) +
                ylab("Inverse simpson index") +
                scale_y_continuous(breaks = seq(0, 30, 10), 
                                   limits = c(0,30), expand = c(0,0))

# richness
# chao plots
ca.chao.p <- ggboxplot(alpha.meta.caeca, x="dpi", y="chao", 
                fill="group2", palette=c("deepskyblue3", "white"), 
              legend.title = "Group") +
              theme(legend.position="none", axis.title.x=element_blank()) +
              ylab("Chao index") + xlab("Days post-treatment") +
              scale_y_continuous(breaks = seq(0, 1000, 200), 
                                 limits = c(0,1000), expand = c(0,0)) +
              stat_compare_means(aes(group = group2), method = "wilcox.test", 
                                 label = "p.signif", hide.ns = TRUE, 
                                    label.y = 900, size=10)




## ileum ##

# diversity
# Inverse simpson plots
il.invsim.p <- ggboxplot(alpha.meta.ileum, x="dpi", y="invsimpson", 
                  fill="group2", palette=c("deepskyblue3", "white")) + 
                theme(legend.position="none") +
                ylab("Inverse simpson index") + xlab("Days post-treatment") +
                scale_y_continuous(breaks = seq(0, 30, 10), 
                                   limits = c(0,30), expand = c(0,0))
                #stat_compare_means(aes(group = chem_administration), method = "t.test", 
                  #label = "p.signif", label.y = 15, size=4)

# richness
# chao plots
il.chao.p <- ggboxplot(alpha.meta.ileum, x="dpi", y="chao", 
                fill="group2", palette=c("deepskyblue3", "white"), 
              legend.title = "Group") +
              theme(legend.position="right") +
              ylab("Chao index") + xlab("Days post-treatment") +
              scale_y_continuous(breaks = seq(0, 1000, 200), 
                                 limits = c(0, 1000), expand = c(0,0)) 





# print .tiff
tiff(filename="results/figures/fig2_alpha.tiff", width=250, 
      height=200, units="mm", res=300)
ggdraw() +
  draw_plot(ca.invsim.p, x=0.03, y=0.5, width=0.41, height=0.45) +
  draw_plot(ca.chao.p, x=0.49, y=0.5, width=0.42, height=0.45) +
  draw_plot(il.invsim.p, x=0.03, y=0, width=0.41, height=0.45) +
  draw_plot(il.chao.p, x=0.49, y=0, width=0.51, height=0.45) +
  draw_plot_label(label = c("A", "B", "C", "D"), size=15, 
                  x = c(0, 0.48, 0, 0.48), y = c(1, 1, 0.5, 0.5))
dev.off()



