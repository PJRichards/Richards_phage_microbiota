######################################
#
# Figure S7
# Histogram of LEfSE results for
# Campylobacter-colonized birds
#
######################################

library("ggplot2")
library("tidyr")
library("cowplot")

# read in LEFSE data
lefse.ca.1dpi <- as.matrix(read.table("data/mothur/phage_ecol.caec_d1.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))
#lefse.ca.2dpi <- as.matrix(read.table("data/mothur/phage_ecol.caec_d2.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE)) # No responsive OTU
lefse.ca.3dpi <- as.matrix(read.table("data/mothur/phage_ecol.caec_d3.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))
lefse.ca.4dpi <- as.matrix(read.table("data/mothur/phage_ecol.caec_d4.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))
lefse.ca.5dpi <- as.matrix(read.table("data/mothur/phage_ecol.caec_d5.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))

# read in taxonomy
tax_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy", header=TRUE)
tax_table <- separate(tax_table, Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), ";")

# this obscures 'rare.otu's'
lefse.ca.1dpi.tax <- merge(lefse.ca.1dpi, tax_table, by="OTU")
#lefse.ca.2dpi.tax <- merge(lefse.ca.2dpi, tax_table, by="OTU")
lefse.ca.3dpi.tax <- merge(lefse.ca.3dpi, tax_table, by="OTU")
lefse.ca.4dpi.tax <- merge(lefse.ca.4dpi, tax_table, by="OTU")
lefse.ca.5dpi.tax <- merge(lefse.ca.5dpi, tax_table, by="OTU")

lefse.ca.1dpi.tax.clean <- lefse.ca.1dpi.tax[complete.cases(lefse.ca.1dpi.tax), ]
lefse.ca.1dpi.tax.clean$ID <- paste(lefse.ca.1dpi.tax.clean$OTU, "..", lefse.ca.1dpi.tax.clean$genus)
lefse.ca.1dpi.tax.clean$fill <- ifelse(lefse.ca.1dpi.tax.clean$Class=="cmp", "deepskyblue4", 
                                  ifelse(lefse.ca.1dpi.tax.clean$Class == "cmpphg", "white", NA))

#lefse.ca.2dpi.tax.clean <- lefse.ca.2dpi.tax[complete.cases(lefse.ca.2dpi.tax), ]
#lefse.ca.2dpi.tax.clean$ID <- paste(lefse.ca.2dpi.tax.clean$OTU, "..", lefse.ca.2dpi.tax.clean$genus)

lefse.ca.3dpi.tax.clean <- lefse.ca.3dpi.tax[complete.cases(lefse.ca.3dpi.tax), ]
lefse.ca.3dpi.tax.clean$ID <- paste(lefse.ca.3dpi.tax.clean$OTU, "..", lefse.ca.3dpi.tax.clean$genus)
lefse.ca.3dpi.tax.clean$fill <- ifelse(lefse.ca.3dpi.tax.clean$Class=="cmp", "deepskyblue4", 
                                  ifelse(lefse.ca.3dpi.tax.clean$Class == "cmpphg", "white", NA))


lefse.ca.4dpi.tax.clean <- lefse.ca.4dpi.tax[complete.cases(lefse.ca.4dpi.tax), ]
lefse.ca.4dpi.tax.clean$ID <- paste(lefse.ca.4dpi.tax.clean$OTU, "..", lefse.ca.4dpi.tax.clean$genus)
lefse.ca.4dpi.tax.clean$fill <- ifelse(lefse.ca.4dpi.tax.clean$Class=="cmp", "deepskyblue4", 
                                  ifelse(lefse.ca.4dpi.tax.clean$Class == "cmpphg", "white", NA))


lefse.ca.5dpi.tax.clean <- lefse.ca.5dpi.tax[complete.cases(lefse.ca.5dpi.tax), ]
lefse.ca.5dpi.tax.clean$ID <- paste(lefse.ca.5dpi.tax.clean$OTU, "..", lefse.ca.5dpi.tax.clean$genus)
lefse.ca.5dpi.tax.clean$fill <- ifelse(lefse.ca.5dpi.tax.clean$Class=="cmp", "deepskyblue4", 
                                  ifelse(lefse.ca.5dpi.tax.clean$Class == "cmpphg", "white", NA))



lefse.ca.1dpi.p <- ggplot(lefse.ca.1dpi.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), fill=fill)) +
                    geom_bar(stat="identity", position="dodge", colour="black") + 
                    scale_fill_identity(guide="legend", labels="Group Cj_phg") +
                    theme(panel.background=element_blank(), panel.grid.major.x=element_line("lightgrey", size = 0.25),
                      legend.text=element_text(size=13),legend.position="right", 
                      legend.title=element_blank(), axis.line=element_line(colour="black", size=0.7)) +
                    ylab("LDA") + xlab("OTU / genus") +
                    scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), expand = c(0,0)) +
                    coord_flip()

#lefse.ca.2dpi.p <- ggplot(lefse.ca.2dpi.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), fill=Class)) +
#  geom_bar(stat="identity", position="dodge", colour="black") + 
#  scale_fill_manual(values=c("white", "deepskyblue4")) +
#  theme(panel.background=element_blank(), panel.grid.major.x=element_line("lightgrey", size = 0.25),
#        legend.text=element_text(size=13),legend.position="none", 
#        legend.title=element_blank(), axis.line=element_line(colour="black", size=0.7)) +
#  ylab("LDA") + xlab("OTU / genus") +
#  scale_y_continuous(breaks = seq(0, 4.5, 1), limits = c(0,4.5), expand = c(0,0)) +
#  coord_flip()

lefse.ca.3dpi.p <- ggplot(lefse.ca.3dpi.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), fill=fill)) +
                    geom_bar(stat="identity", position="dodge", colour="black") + 
                    scale_fill_identity(guide="legend", labels="Group Cj") +
                    theme(panel.background=element_blank(), panel.grid.major.x=element_line("lightgrey", size = 0.25),
                      legend.text=element_text(size=13),legend.position="right", 
                      legend.title=element_blank(), axis.line=element_line(colour="black", size=0.7)) +
                    ylab("LDA") + xlab("OTU / genus") +
                    scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), expand = c(0,0)) +
                    coord_flip()

lefse.ca.4dpi.p <- ggplot(lefse.ca.4dpi.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), fill=fill)) +
                    geom_bar(stat="identity", position="dodge", colour="black") + 
                    scale_fill_identity(guide="legend", labels="Group Cj") +
                    theme(panel.background=element_blank(), panel.grid.major.x=element_line("lightgrey", size = 0.25),
                      legend.text=element_text(size=13),legend.position="right", 
                      legend.title=element_blank(), axis.line=element_line(colour="black", size=0.7)) +
                    ylab("LDA") + xlab("OTU / genus") +
                    scale_y_continuous(breaks = seq(0, 4.5, 1), limits = c(0,4.5), expand = c(0,0)) +
                    coord_flip()

lefse.ca.5dpi.p <- ggplot(lefse.ca.5dpi.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), fill=fill)) +
                    geom_bar(stat="identity", position="dodge", colour="black") + 
                    scale_fill_identity(guide="legend", labels="Group Cj_phg") +
                    theme(panel.background=element_blank(), panel.grid.major.x=element_line("lightgrey", size = 0.25),
                      legend.text=element_text(size=13), legend.position="right",
                      legend.title=element_blank(), axis.line=element_line(colour="black", size=0.7)) +
                    ylab("LDA") + xlab("OTU / genus") +
                    scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), expand = c(0,0)) +
                    coord_flip()

tiff(filename="results/figures/figS7_ca_LEfSE.tiff", width=250, height=250, units="mm", res=300)
ggdraw() +
  draw_plot(lefse.ca.1dpi.p, x=0.05, y=0.78, width=0.8, height=0.22) +
  draw_plot(lefse.ca.3dpi.p, x=0.05, y=0.55, width=0.8, height=0.22) +
  draw_plot(lefse.ca.4dpi.p, x=0.05, y=0.32, width=0.8, height=0.22) +
  draw_plot(lefse.ca.5dpi.p, x=0.05, y=0, width=0.8, height=0.31) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20, x = c(0, 0, 0, 0), y = c(0.99, 0.77, 0.53, 0.30))
dev.off()
