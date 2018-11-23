####################################################
#
# Fig. S5. 
# Phyla-level stacked barcharts
#
####################################################

library("ggplot2")
library("data.table")
library("reshape")
library("tidyr")
library("cowplot")

# read in and format data
OTU_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared", header=TRUE, row.names=2)
OTU_table <- OTU_table[-c(1,2)] # tidy

shallow.community.remove <- c("il_cmp_d2_2", "il_cmpphg_d4_3")
OTU_table <- OTU_table[!(row.names(OTU_table) %in% shallow.community.remove), ]


tax_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy", header=TRUE, row.names=1)
tax_table <- separate(tax_table, Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), ";")

# subset by site
OTU_table.ca <- OTU_table[grep("^ca", row.names(OTU_table)),]
OTU_table.il <- OTU_table[grep("^il", row.names(OTU_table)),]


# remove OTUs reported in only 2 communities
# appropriate for phyla overview? 
OTU_table.ca.clean <- OTU_table.ca[apply(OTU_table.ca, 2, function(x)sum(x>0)) > 1]
OTU_table.il.clean <- OTU_table.il[apply(OTU_table.il, 2, function(x)sum(x>0)) > 1]


# calculate relative abundance 
OTU_table.ca.RA <- 100*(OTU_table.ca.clean/rowSums(OTU_table.ca.clean))
OTU_table.il.RA <- 100*(OTU_table.il.clean/rowSums(OTU_table.il.clean))



# transpose OTUs
OTU_table.ca.RA.t <- as.data.frame(t(OTU_table.ca.RA))
OTU_table.il.RA.t <- as.data.frame(t(OTU_table.il.RA))



# merge taxonomy and OTU
OTU_table.ca.RA.tax <- merge(OTU_table.ca.RA.t, tax_table, by=0) # merge by row.name
OTU_table.ca.RA.tax$Row.names <- NULL

OTU_table.il.RA.tax <- merge(OTU_table.il.RA.t, tax_table, by=0)
OTU_table.il.RA.tax$Row.names <- NULL




# summarize as phyla
OTU_table.ca.RA.phyla <- aggregate(OTU_table.ca.RA.tax[, 1:52], list(OTU_table.ca.RA.tax$phylum), sum)
OTU_table.il.RA.phyla <- aggregate(OTU_table.il.RA.tax[, 1:42], list(OTU_table.il.RA.tax$phylum), sum)




# reformat as data.frame
OTU_table.ca.RA.phyla.t <- as.data.frame(t(OTU_table.ca.RA.phyla))
colnames(OTU_table.ca.RA.phyla.t) <- c("Actinobacteria(100)", "Bacteria_unclassified(100)", 
                                       "Bacteroidetes(100)", "Firmicutes(100)", "Proteobacteria(100)")
OTU_table.ca.RA.phyla.t <- OTU_table.ca.RA.phyla.t[-1,]
setDT(OTU_table.ca.RA.phyla.t, keep.rownames = TRUE)[]
colnames(OTU_table.ca.RA.phyla.t)[1] <- "group"



OTU_table.il.RA.phyla.t <- as.data.frame(t(OTU_table.il.RA.phyla))
colnames(OTU_table.il.RA.phyla.t) <- c("Acidobacteria(100)", "Actinobacteria(100)", "Bacteria_unclassified(100)",
                                       "Bacteroidetes(100)", "Chloroflexi(100)", "Firmicutes(100)", 
                                       "Ignavibacteriae(100)", "Planctomycetes(100)", "Proteobacteria(100)",
                                       "Verrucomicrobia(100)")
OTU_table.il.RA.phyla.t <- OTU_table.il.RA.phyla.t[-1,]
setDT(OTU_table.il.RA.phyla.t, keep.rownames = TRUE)[]
colnames(OTU_table.il.RA.phyla.t)[1] <- "group"

# add in grouping info
OTU_table.ca.RA.phyla.t <- separate(OTU_table.ca.RA.phyla.t, group, 
                                    c("site", "treatment", "time", "replicate"), "_")
OTU_table.ca.RA.phyla.t$ID <- paste(OTU_table.ca.RA.phyla.t$treatment, OTU_table.ca.RA.phyla.t$time, sep = "_")

OTU_table.il.RA.phyla.t <- separate(OTU_table.il.RA.phyla.t, group, 
                                    c("site", "treatment", "time", "replicate"), "_")
OTU_table.il.RA.phyla.t$ID <- paste(OTU_table.il.RA.phyla.t$treatment, OTU_table.il.RA.phyla.t$time, sep = "_")


# pre-melt format
OTU_table.ca.RA.phyla.t$site <- NULL
OTU_table.ca.RA.phyla.t$time <- NULL
OTU_table.ca.RA.phyla.t$replicate <- NULL
OTU_table.ca.RA.phyla.t$treatment <- NULL

OTU_table.ca.RA.phyla.t$`Actinobacteria(100)` <- as.numeric(as.character(OTU_table.ca.RA.phyla.t$`Actinobacteria(100)`))
OTU_table.ca.RA.phyla.t$`Bacteria_unclassified(100)` <- as.numeric(as.character(OTU_table.ca.RA.phyla.t$`Bacteria_unclassified(100)`))
OTU_table.ca.RA.phyla.t$`Bacteroidetes(100)` <- as.numeric(as.character(OTU_table.ca.RA.phyla.t$`Bacteroidetes(100)`))
OTU_table.ca.RA.phyla.t$`Firmicutes(100)` <- as.numeric(as.character(OTU_table.ca.RA.phyla.t$`Firmicutes(100)`))
OTU_table.ca.RA.phyla.t$`Proteobacteria(100)` <- as.numeric(as.character(OTU_table.ca.RA.phyla.t$`Proteobacteria(100)`))

OTU_table.il.RA.phyla.t$site <- NULL
OTU_table.il.RA.phyla.t$time <- NULL
OTU_table.il.RA.phyla.t$replicate <- NULL
OTU_table.il.RA.phyla.t$treatment <- NULL

OTU_table.il.RA.phyla.t$`Acidobacteria(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Acidobacteria(100)`))
OTU_table.il.RA.phyla.t$`Actinobacteria(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Actinobacteria(100)`))
OTU_table.il.RA.phyla.t$`Bacteria_unclassified(100)` <- as.numeric(as.character(
                                                            OTU_table.il.RA.phyla.t$`Bacteria_unclassified(100)`))
OTU_table.il.RA.phyla.t$`Bacteroidetes(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Bacteroidetes(100)`))
OTU_table.il.RA.phyla.t$`Chloroflexi(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Chloroflexi(100)`))
OTU_table.il.RA.phyla.t$`Firmicutes(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Firmicutes(100)`))
OTU_table.il.RA.phyla.t$`Ignavibacteriae(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Ignavibacteriae(100)`))
OTU_table.il.RA.phyla.t$`Planctomycetes(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Planctomycetes(100)`))
OTU_table.il.RA.phyla.t$`Proteobacteria(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Proteobacteria(100)`))
OTU_table.il.RA.phyla.t$`Verrucomicrobia(100)` <- as.numeric(as.character(OTU_table.il.RA.phyla.t$`Verrucomicrobia(100)`))


# summarize by cohort
OTU_table.ca.RA.phyla.mean <- aggregate(OTU_table.ca.RA.phyla.t[, 1:5], list(OTU_table.ca.RA.phyla.t$ID), mean)
OTU_table.il.RA.phyla.mean <- aggregate(OTU_table.il.RA.phyla.t[, 1:10], list(OTU_table.il.RA.phyla.t$ID), mean)

# remove controls
OTU_table.ca.RA.phyla.mean.exp <- OTU_table.ca.RA.phyla.mean[-11,]
OTU_table.il.RA.phyla.mean.exp <- OTU_table.il.RA.phyla.mean[-9,]


# melt data
OTU_table.ca.RA.phyla.mean.exp.m <- melt.data.frame(OTU_table.ca.RA.phyla.mean.exp, id.vars="Group.1")


OTU_table.il.RA.phyla.mean.exp.m <- melt.data.frame(OTU_table.il.RA.phyla.mean.exp, id.vars="Group.1")


# plot graphing

p.ca.phyla <- ggplot(aes(x=Group.1, y=value, fill=variable), data=OTU_table.ca.RA.phyla.mean.exp.m) +
                geom_bar(stat="identity", colour="white") +
                scale_fill_manual(values=c("#FFCC99", "#808080", "#2BCE48", "#0075DC", "#C20088"), 
                  breaks=c("Actinobacteria(100)", "Bacteria_unclassified(100)", "Bacteroidetes(100)", 
                  "Firmicutes(100)", "Proteobacteria(100)")) +
                scale_x_discrete(labels=c("cmp_d1" = "Cj 1 dpt", "cmpphg_d1" = "Cj_phg 1 dpt", "cmp_d2" = "Cj 2 dpt", 
                  "cmpphg_d2" = "Cj_phg 2 dpt", "cmp_d3" = "Cj 3 dpt", "cmpphg_d3" = "Cj_phg 3 dpt", "cmp_d4" = "Cj 4 dpt",
                  "cmpphg_d4" = "Cj_phg 4 dpt", "cmp_d5" = "Cj 5 dpt", "cmpphg_d5" = "Cj_phg 5 dpt")) +
                theme(panel.background=element_blank(), panel.grid.major = element_blank(), 
                  axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.title.y=element_blank(),
                  legend.position="none", legend.title=element_blank(), 
                  axis.line = element_line(colour = 'black', size=0.7)) +
                xlab("Sample")

p.il.phyla <- ggplot(aes(x=Group.1, y=value, fill=variable), data=OTU_table.il.RA.phyla.mean.exp.m) +
                geom_bar(stat="identity", colour="white") +
                scale_fill_manual(values=c("#993F00", "#FFCC99", "#808080", "#2BCE48", "#4C005C", "#0075DC", "#94FFB5",
                  "#9DCC00", "#C20088", "#003380"), breaks=c("Acidobacteria(100)", "Actinobacteria(100)", 
                  "Bacteria_unclassified(100)", "Bacteroidetes(100)", "Chloroflexi(100)", "Firmicutes(100)", 
                  "Ignavibacteriae(100)", "Planctomycetes(100)", "Proteobacteria(100)", "Verrucomicrobia(100)")) +
                scale_x_discrete(labels=c("cmp_d2" = "Cj 2 dpt", "cmpphg_d2" = "Cj_phg 2 dpt", "cmp_d3" = "Cj 3 dpt", 
                  "cmpphg_d3" = "Cj_phg 3 dpt", "cmp_d4" = "Cj 4 dpt", "cmpphg_d4" = "Cj_phg 4 dpt", 
                  "cmp_d5" = "Cj 5 dpt", "cmpphg_d5" = "Cj_phg 5 dpt")) +
                theme(panel.background=element_blank(), panel.grid.major = element_blank(), 
                  axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                  legend.position="left", legend.title=element_blank(), 
                  axis.line = element_line(colour = 'black', size=0.7)) +
                ylab("Relative abundance (%)") +
                xlab("Sample")

# print plots
tiff(filename="results/figures/figS5_phyla.tiff", width=200, height=110, units="mm", res=300)
ggdraw() +
  draw_plot(p.il.phyla, x=0.01, y=0, width=0.63, height=0.97) +
  draw_plot(p.ca.phyla, x=0.65, y=0, width=0.34, height=0.97) +
  draw_plot_label(label = c("A", "B"), size=12, x=c(0.365, 0.64), y=c(1, 1))
dev.off()
  