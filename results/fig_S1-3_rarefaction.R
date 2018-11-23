###############################################
#
# Figure S1.
# Rarefaction curves
#
##############################################

# read in data
rarefaction <- read.table(file="data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.groups.rarefaction", sep = "\t", header=TRUE)



### Ceca ###
# rarefaction curves for caecal communities from group Cj and group Cj_phg cohorts.


png(filename="results/figures/figS1_exp_caeca_rarefaction.tiff", width=100, height=160, units="mm", res=300)
par(mfrow=c(5,2))


# 1 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:672], y=rarefaction$X0.03.ca_cmp_d1_1[1:672], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,400), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d1_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d1_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d1_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d1_5, type="l", lwd=0.2) # maybe fail?
title(main="A i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:672], y=rarefaction$X0.03.ca_cmpphg_d1_1[1:672], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d1_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d1_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d1_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d1_5, type="l", lwd=0.2)
title(main="A ii.", adj=0, line=0.5, cex.main=0.9)



# 2 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:1402], y=rarefaction$X0.03.ca_cmp_d2_1[1:1402], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d2_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d2_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d2_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d2_5, type="l", lwd=0.2)
title(main="B i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:1402], y=rarefaction$X0.03.ca_cmpphg_d2_1[1:1402], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,600), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d2_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d2_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d2_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d2_5, type="l", lwd=0.2)
title(main="B ii.", adj=0, line=0.5, cex.main=0.9)


# 3 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:1865], y=rarefaction$X0.03.ca_cmp_d3_1[1:1865], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,600), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d3_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d3_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d3_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d3_5, type="l", lwd=0.2)
title(main="C i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:1865], y=rarefaction$X0.03.ca_cmpphg_d3_1[1:1865], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,600), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d3_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d3_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d3_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d3_5, type="l", lwd=0.2)
title(main="C ii.", adj=0, line=0.5, cex.main=0.9)


# 4 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:1865], y=rarefaction$X0.03.ca_cmp_d4_1[1:1865], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d4_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d4_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d4_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d4_5, type="l", lwd=0.2)
title(main="D i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:1865], y=rarefaction$X0.03.ca_cmpphg_d4_1[1:1865], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d4_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d4_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d4_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d4_5, type="l", lwd=0.2)
title(main="D ii.", adj=0, line=0.5, cex.main=0.9)


# 5 dpi #

# cmp
# NB. no bird 5
par(mar=c(4, 5, 1.1, 0.5))
plot(x=rarefaction$numsampled[1:1098], y=rarefaction$X0.03.ca_cmp_d5_1[1:1098], 
      xlab="Number of reads sampled", ylab="OTUs", type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d5_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d5_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmp_d5_4, type="l", lwd=0.2) 
title(main="E i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(4, 2.5, 1.1, 3))
plot(x=rarefaction$numsampled[1:1098], y=rarefaction$X0.03.ca_cmpphg_d5_1[1:1098], 
      xlab="Number of reads sampled", ylab="OTUs", type="l", ylim=c(0,500), lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d5_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d5_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d5_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_cmpphg_d5_5, type="l", lwd=0.2)
title(main="E ii.", adj=0, line=0.5, cex.main=0.9)


dev.off()



### Controls ###
# rarefaction curves for caecal communities from ctl cohort.

png(filename="results/figures/figS2_ctl_caeca_rarefaction.tiff", width=100, height=160, units="mm", res=300)
par(mfrow=c(5,2))

par(mar=c(4, 5, 1.1, 0.5))
plot(x=rarefaction$numsampled[1:1098], y=rarefaction$X0.03.ca_ctl_d5_1[1:1098], 
     xlab="Number of reads sampled", ylab=NA, type="l", ylim=c(0,400), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_ctl_d5_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.ca_ctl_d5_3, type="l", lwd=0.2)
title(adj=0, line=0.5, cex.main=0.9)


dev.off()



### Ileum ###
# rarefaction curves for ileal communities from group Cj and group Cj_phg cohorts.

png(filename="results/figures/figS3_ileum_rarefaction.tiff", width=100, height=160, units="mm", res=300)
par(mfrow=c(5,2))


# 2 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:4393], y=rarefaction$X0.03.il_cmp_d2_1[1:4393], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,200), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d2_2, type="l", lwd=0.2) # **fail**
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d2_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d2_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d2_5, type="l", lwd=0.2)
title(main="A i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:4393], y=rarefaction$X0.03.il_cmpphg_d2_1[1:4393], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,200), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d2_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d2_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d2_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d2_5, type="l", lwd=0.2)
title(main="A ii.", adj=0, line=0.5, cex.main=0.9)


# 3 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:2616], y=rarefaction$X0.03.il_cmp_d3_1[1:2616], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,400), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d3_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d3_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d3_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d3_5, type="l", lwd=0.2)
title(main="B i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:2616], y=rarefaction$X0.03.il_cmpphg_d3_1[1:2616], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,400), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d3_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d3_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d3_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d3_5, type="l", lwd=0.2)
title(main="B ii.", adj=0, line=0.5, cex.main=0.9)


# 4 dpi #

# cmp
par(mar=c(3.3, 5, 1.4, 0.5))
plot(x=rarefaction$numsampled[1:1098], y=rarefaction$X0.03.il_cmp_d4_1[1:1098], 
     xlab=NA, ylab="OTUs", type="l", ylim=c(0,100), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d4_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d4_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d4_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d4_5, type="l", lwd=0.2)
title(main="C i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(3.3, 2.5, 1.4, 3))
plot(x=rarefaction$numsampled[1:4117], y=rarefaction$X0.03.il_cmpphg_d4_1[1:4117], 
     xlab=NA, ylab=NA, type="l", ylim=c(0,900), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d4_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d4_3, type="l", lwd=0.2) #**fail**
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d4_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d4_5, type="l", lwd=0.2)
title(main="C ii.", adj=0, line=0.5, cex.main=0.9)


# 5 dpi #

# cmp
par(mar=c(4, 5, 1.1, 0.5))
plot(x=rarefaction$numsampled[1:2616], y=rarefaction$X0.03.il_cmp_d5_1[1:2616], 
     xlab="Number of reads sampled", ylab="OTUs", type="l", ylim=c(0,100), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d5_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d5_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d5_4, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmp_d5_5, type="l", lwd=0.2)
title(main="D i.", adj=0, line=0.5, cex.main=0.9)

# cmpphg
par(mar=c(4, 2.5, 1.1, 3))
plot(x=rarefaction$numsampled[1:2616], y=rarefaction$X0.03.il_cmpphg_d5_1[1:2616], 
     xlab="Number of reads sampled", ylab="OTUs", type="l", ylim=c(0,200), lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d5_2, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d5_3, type="l", lwd=0.2)
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d5_4, type="l", lwd=0.2) 
points(x=rarefaction$numsampled, y=rarefaction$X0.03.il_cmpphg_d5_5, type="l", lwd=0.2)
title(main="D ii.", adj=0, line=0.5, cex.main=0.9)

dev.off()
