############################################################
#
# Figure S8. LEfSE analysis for OTUs discriminatory of 
# *Campylobacter*-colonized and *Campylobacter*-free birds
#
############################################################

library("ggplot2")
library("tidyr")
library("cowplot")

## read in LEFSE data

# ctl vs cmp vs cmpphg 3-way comparison
lefse.ca.all <- as.matrix(read.table("data/mothur/phage_ecol.caec_d5_all.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))

# ctl vs cmp
lefse.ca.cmp_v_ctl <- as.matrix(read.table("data/mothur/phage_ecol.caec_d5_cmp.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))

# ctl vs cmpphg
lefse.ca.cmpphg_v_ctl <- as.matrix(read.table("data/mothur/phage_ecol.caec_d5_cmpphg.0.03.filter.0.03.subsample.0.03.lefse_summary", fill=T, header = TRUE))

## read in taxonomy
tax_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy", header=TRUE)
tax_table <- separate(tax_table, Taxonomy, 
              c("kingdom", "phylum", "class", "order", "family", "genus"), ";")


## obscure 'rare.otu's'


# all
lefse.ca.all.tax <- merge(lefse.ca.all, tax_table, by="OTU")
lefse.ca.all.tax.clean <- lefse.ca.all.tax[complete.cases(lefse.ca.all.tax), ]
lefse.ca.all.tax.clean$ID <- paste(lefse.ca.all.tax.clean$OTU, "..", lefse.ca.all.tax.clean$genus)
lefse.ca.all.tax.clean$fill <- ifelse(lefse.ca.all.tax.clean$Class=="cmp", "deepskyblue4",
                                      ifelse(lefse.ca.all.tax.clean$Class=="cmpphg", "white", 
                                        NA))

# cmp_v_ctl
lefse.ca.cmp_v_ctl.tax <- merge(lefse.ca.cmp_v_ctl, tax_table, by="OTU")
lefse.ca.cmp_v_ctl.tax.clean <- lefse.ca.cmp_v_ctl.tax[complete.cases(lefse.ca.cmp_v_ctl.tax), ]
lefse.ca.cmp_v_ctl.tax.clean$ID <- paste(lefse.ca.cmp_v_ctl.tax.clean$OTU, "..", 
                                      lefse.ca.cmp_v_ctl.tax.clean$genus)
lefse.ca.cmp_v_ctl.tax.clean$fill <- ifelse(lefse.ca.cmp_v_ctl.tax.clean$Class=="cmp", 
                                      "deepskyblue4", ifelse(lefse.ca.cmp_v_ctl.tax.clean$Class
                                        =="ctl", "red3", NA))

# cmpphg_v_ctl
lefse.ca.cmpphg_v_ctl.tax <- merge(lefse.ca.cmpphg_v_ctl, tax_table, by="OTU")
lefse.ca.cmpphg_v_ctl.tax.clean <- lefse.ca.cmpphg_v_ctl.tax[complete.cases(
                                      lefse.ca.cmpphg_v_ctl.tax), ]
lefse.ca.cmpphg_v_ctl.tax.clean$ID <- paste(lefse.ca.cmpphg_v_ctl.tax.clean$OTU, "..", 
                                        lefse.ca.cmpphg_v_ctl.tax.clean$genus)
lefse.ca.cmpphg_v_ctl.tax.clean$fill <- ifelse(lefse.ca.cmpphg_v_ctl.tax.clean$Class=="cmp",
                                                "deepskyblue4", 
                                                  ifelse(lefse.ca.cmpphg_v_ctl.tax.clean$Class
                                                    =="cmpphg", "white", NA))


lefse.ca.all.p <- ggplot(lefse.ca.all.tax.clean, aes(x=ID, y=as.numeric(as.character(LDA)), 
                    fill=fill)) +
                      geom_bar(stat="identity", position="dodge", colour="black") + 
                      scale_fill_identity(guide="legend", labels=c("Group Cj_phg", "Control")) +
                      theme(panel.background=element_blank(), 
                      panel.grid.major.x=element_line("lightgrey", size = 0.25),
                        legend.text=element_text(size=13),legend.position="right", 
                          legend.title=element_blank(), axis.line=element_line(colour="black", 
                            size=0.7)) +
                      ylab("LDA") + xlab("OTU / genus") +
                      scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), expand = c(0,0)) +
                      coord_flip()


lefse.ca.cmp_v_ctl.p <- ggplot(lefse.ca.cmp_v_ctl.tax.clean, aes(x=ID, 
                          y=as.numeric(as.character(LDA)), fill=fill)) +
                            geom_bar(stat="identity", position="dodge", colour="black") + 
                            scale_fill_identity(guide="legend", labels=c("Group Cj", "Control")) +
                            theme(panel.background=element_blank(), 
                              panel.grid.major.x=element_line("lightgrey", size = 0.25),
                                legend.text=element_text(size=13),legend.position="right", 
                                  legend.title=element_blank(), 
                                    axis.line=element_line(colour="black", size=0.7)) +
                            ylab("LDA") + xlab("OTU / genus") +
                            scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), 
                              expand = c(0,0)) +
                            coord_flip()


lefse.ca.cmpphg_v_ctl.p <- ggplot(lefse.ca.cmpphg_v_ctl.tax.clean, aes(x=ID, 
                            y=as.numeric(as.character(LDA)), fill=fill)) +
                              geom_bar(stat="identity", position="dodge", colour="black") + 
                              scale_fill_identity(guide="legend", labels="Group Cj_phg") +
                              theme(panel.background=element_blank(), 
                                panel.grid.major.x=element_line("lightgrey", size = 0.25),
                                legend.text=element_text(size=13),legend.position="right", 
                                legend.title=element_blank(), 
                                axis.line=element_line(colour="black", size=0.7)) +
                              ylab("LDA") + xlab("OTU / genus") +
                              scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), 
                              expand = c(0,0)) +
                              coord_flip()



tiff(filename="results/figures/figS8_ctl_LEfSE.tiff", width=250, height=250, units="mm", res=300)
ggdraw() +
  draw_plot(lefse.ca.all.p, x=0.05, y=0.66, width=0.9, height=0.32) +
  draw_plot(lefse.ca.cmp_v_ctl.p, x=0.05, y=0.33, width=0.9, height=0.32) +
  draw_plot(lefse.ca.cmpphg_v_ctl.p, x=0.05, y=0, width=0.9, height=0.32) +
  draw_plot_label(label = c("A", "B", "C"), size = 20, x = c(0, 0, 0), y = c(0.99, 0.66, 0.33))
dev.off()
