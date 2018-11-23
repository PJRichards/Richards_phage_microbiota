###############################
#
# Figure 1
# Campy CFU. 
#
###############################

library("ggpubr")
library("cowplot")
library("plyr")

# read in data
campy.cfu <- read.table(file="data/references/campy_cfu.txt", sep="\t", header=TRUE, row.names=1)
campy_phage.pfu <- read.table(file="data/references/campy_phage_pfu.txt", sep="\t", header=TRUE)


# tidy CFU
campy.cfu$cfu_ml <- campy.cfu$dil_factor*rowMeans(subset(campy.cfu, select = c(cfu_count_1, cfu_count_2, cfu_count_3)))
campy.cfu$log_cfu_ml <- log10(campy.cfu$cfu_ml)
campy.cfu$dpi <- as.factor(campy.cfu$dpi)
campy.cfu$group <- as.factor(campy.cfu$group)


## tidy pfu
campy_phage.pfu$dil_factor <- revalue(campy_phage.pfu$dil_factor, c("not applicable"="0"))
campy_phage.pfu$dil_factor <- as.numeric(as.character(campy_phage.pfu$dil_factor))
campy_phage.pfu$dpi <- as.factor(campy_phage.pfu$dpi)
campy_phage.pfu$group <- as.factor(campy_phage.pfu$group)

# cp20
campy_phage.pfu$cp20_pfu_g <- campy_phage.pfu$dil_factor*rowMeans(subset(campy_phage.pfu, select = c(cp20_count_1, cp20_count_2, cp20_count_3)))
campy_phage.pfu$cp20_pfu_g <- revalue(as.character(campy_phage.pfu$cp20_pfu_g), c("0"="1000"))
campy_phage.pfu$cp20_pfu_g <- as.numeric(campy_phage.pfu$cp20_pfu_g)
campy_phage.pfu$log_cp20_pfu_g <- log10(campy_phage.pfu$cp20_pfu_g)

# cp30
campy_phage.pfu$cp30_pfu_g <- campy_phage.pfu$dil_factor*rowMeans(subset(campy_phage.pfu, select = c(cp30_count_1, cp30_count_2, cp30_count_3)))
campy_phage.pfu$cp30_pfu_g <- revalue(as.character(campy_phage.pfu$cp30_pfu_g), c("0"="1000"))
campy_phage.pfu$cp30_pfu_g <- as.numeric(campy_phage.pfu$cp30_pfu_g)
campy_phage.pfu$log_cp30_pfu_g <- log10(campy_phage.pfu$cp30_pfu_g)



# subset
campy.cfu.caeca <- campy.cfu[campy.cfu$site=="caeca",]
campy.cfu.ileum <- campy.cfu[campy.cfu$site=="ileum",]
campy.cfu.colon <- campy.cfu[campy.cfu$site=="colon",]

campy.pfu.caeca <- campy_phage.pfu[campy_phage.pfu$site=="caeca",]
campy.pfu.ileum <- campy_phage.pfu[campy_phage.pfu$site=="ileum",]



ca.cfu.p <- ggdotplot(campy.cfu.caeca, x="dpi", y="log_cfu_ml", fill="group", size=0.5) +
  xlab("Days post-treatment") +
  ylab(expression(italic("Campylobacter")~textstyle("counts (log"[10]~textstyle("CFU g"^-1*textstyle(")"))))) +
          theme(legend.position=c(0.9, 0.15), legend.box.background = element_rect(colour = "black"), 
                legend.background = element_blank()) + 
  scale_fill_manual(values=c("deepskyblue3", "white"), breaks=c("1", "2"), labels=c("Cj", "Cj_phg")) +
  stat_summary(fun.y=mean, geom="point", shape = 15,  aes(group = group), 
          position=position_dodge(0.8), color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", aes(group = group), 
          position=position_dodge(0.8), width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 9, 1), limits = c(2.4,9.4), expand = c(0,0)) +
  stat_compare_means(aes(group = group), method = "t.test", label = "p.signif", 
          hide.ns = TRUE, label.y = 9.2, size=6)

il.cfu.p <- ggdotplot(campy.cfu.ileum, x = "dpi", y = "log_cfu_ml", fill = "group", palette=c("deepskyblue3", "white"), size = 0.8) +
  xlab("Days post-treatment") +
  ylab(expression(atop(italic("Campylobacter")~textstyle("counts"), paste(textstyle("(log"[10]~textstyle("CFU g"^-1*textstyle(")"))))))) +
  theme(legend.position="none", text=element_text(size=8)) + 
  stat_summary(fun.y=mean, geom="point", shape = 15, aes(group = group), 
               position=position_dodge(0.8), color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", aes(group = group), 
               position=position_dodge(0.8), width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 9, 2), limits = c(2.4,9.4), expand = c(0,0)) +
  stat_compare_means(aes(group = group), method = "t.test", label = "p.signif", 
                     hide.ns = TRUE, label.y = 9, size=3)

col.cfu.p <- ggdotplot(campy.cfu.colon, x = "dpi", y = "log_cfu_ml", fill = "group", 
                       palette=c("deepskyblue3", "white"), size = 0.8) +
  xlab("Days post-treatment") +
  theme(legend.position="none", axis.title.y=element_blank(), text=element_text(size=8)) + 
  stat_summary(fun.y=mean, geom="point", shape = 15,  aes(group = group), 
               position=position_dodge(0.8), color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", aes(group = group), 
               position=position_dodge(0.8), width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 9, 2), limits = c(2.4,9.4), expand = c(0,0)) +
  stat_compare_means(aes(group = group), method = "t.test", label = "p.signif", 
                     hide.ns = TRUE, label.y = 9, size=3)

ca.pfu.cp20.p <- ggdotplot(campy.pfu.caeca, x = "dpi", y = "log_cp20_pfu_g ", palette=c("deepskyblue3"), size = 0.8) +
  theme(text=element_text(size=8)) +
  xlab("Days post-treatment") + 
  ylab(expression(atop(textstyle("Phage counts"), paste(textstyle("(log"[10]~textstyle("PFU g"^-1*textstyle(")"))))))) +
  stat_summary(fun.y=mean, geom="point", shape = 15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 9, 2), limits = c(2.4, 8.4), expand = c(0,0))

ca.pfu.cp30.p <- ggdotplot(campy.pfu.caeca, x = "dpi", y = "log_cp30_pfu_g ", palette=c("deepskyblue3"), size = 0.8) +
  theme(axis.title.y=element_blank(), text=element_text(size=8)) +
  xlab("Days post-treatment") + 
  stat_summary(fun.y=mean, geom="point", shape = 15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 9, 2), limits = c(2.4, 8.4), expand = c(0,0))





tiff(filename="results/figures/fig1_counts.tiff", width=190, height=110, units="mm", res=300)
ggdraw() +
  draw_plot(ca.cfu.p, x=0, y=0, width=0.55, height=0.95) +
  draw_plot(il.cfu.p, x=0.56, y=0.01, width=0.25, height=0.4) +
  draw_plot(col.cfu.p, x=0.81, y=0.01, width=0.19, height=0.4) +
  draw_plot(ca.pfu.cp20.p, x=0.56, y=0.53, width=0.25, height=0.4) + 
  draw_plot(ca.pfu.cp30.p, x=0.81, y=0.53, width=0.19, height=0.4) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 12, x = c(0, 0.61, 0.8, 0.61, 0.8), y = c(1, 1, 1, 0.48, 0.48))
dev.off()
  

 


