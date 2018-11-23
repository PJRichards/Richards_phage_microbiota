###############################
#
# Figure S4  
# Phage PFU. 
#
###############################

library("ggpubr")
library("cowplot")
library("plyr")

# read in data
campy_phage.pfu <- read.table(file="data/references/campy_phage_pfu.txt", sep="\t", header=TRUE)


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
campy.pfu.colon <- campy_phage.pfu[campy_phage.pfu$site=="colon",]
campy.pfu.ileum <- campy_phage.pfu[campy_phage.pfu$site=="ileum",]


# make plots
il.pfu.cp20.p <- ggdotplot(campy.pfu.ileum, x="dpi", y="log_cp20_pfu_g ", palette=c("deepskyblue3")) +
  xlab("Days post-treatment") + 
  ylab(expression("CP20 phage counts (log"[10]~textstyle("PFU g"^-1*textstyle(")")))) +
  theme(text=element_text(size=8)) +
  stat_summary(fun.y=mean, geom="point", shape=15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(2.4, 7.4), expand = c(0,0))

il.pfu.cp30.p <- ggdotplot(campy.pfu.ileum, x="dpi", y="log_cp30_pfu_g ", palette=c("deepskyblue3")) +
  xlab("Days post-treatment") + 
  ylab(expression("CP30 phage counts (log"[10]~textstyle("PFU g"^-1*textstyle(")")))) +
  theme(text=element_text(size=8)) +
  stat_summary(fun.y=mean, geom="point", shape=15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(2.4, 7.4), expand = c(0,0))

colon.pfu.cp20.p <- ggdotplot(campy.pfu.colon, x="dpi", y="log_cp20_pfu_g ", palette=c("deepskyblue3")) +
  xlab("Days post-treatment") + 
  ylab(expression("CP20 phage counts (log"[10]~textstyle("PFU g"^-1*textstyle(")")))) +
  theme(text=element_text(size=8)) +
  stat_summary(fun.y=mean, geom="point", shape = 15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(2.4, 7.4), expand = c(0,0))


colon.pfu.cp30.p <- ggdotplot(campy.pfu.colon, x="dpi", y="log_cp30_pfu_g ", palette=c("deepskyblue3")) +
  xlab("Days post-treatment") + 
  ylab(expression("CP30 phage counts (log"[10]~textstyle("PFU g"^-1*textstyle(")")))) +
  theme(text=element_text(size=8)) +
  stat_summary(fun.y=mean, geom="point", shape = 15,  color="red") +
  stat_summary(fun.data=mean_sd, geom="errorbar", width=0.2, color="red") +
  scale_y_continuous(breaks = seq(3, 7, 1), limits = c(2.4, 7.4), expand = c(0,0))



tiff(filename="results/figures/figS4_pfu.tiff", width=125, height=125, units="mm", res=300)
ggdraw() +
  draw_plot(il.pfu.cp20.p, x=0.05, y=0.5, width=0.4, height=0.4) +
  draw_plot(il.pfu.cp30.p, x=0.55, y=0.5, width=0.4, height=0.4) +
  draw_plot(colon.pfu.cp20.p, x=0.05, y=0, width=0.4, height=0.4) +
  draw_plot(colon.pfu.cp30.p, x=0.55, y=0, width=0.4, height=0.4) + 
  draw_plot_label(label = c("A", "B", "C", "D"), size=10, 
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
dev.off()
  

 


