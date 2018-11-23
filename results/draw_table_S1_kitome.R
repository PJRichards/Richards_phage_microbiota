####################################################
#
# Table S1. 
# No template controls
#
####################################################


library("data.table")
library("reshape")
library("tidyr")
library("gridExtra")

# read in and format data
OTU_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared", header = TRUE, row.names = 2)
OTU_table <- OTU_table[-c(1,2)] # tidy

tax_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy", header=TRUE, row.names=1)
tax_table <- separate(tax_table, Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), ";")

# subset by site
OTU_table.ctl <- OTU_table[c(99:104),]

# remove OTUs < 1 read
OTU_table.ctl.clean <- OTU_table.ctl[colSums(OTU_table.ctl) > 1]

# calculate relative abundance 
OTU_table.ctl.RA <- 100*(OTU_table.ctl.clean/rowSums(OTU_table.ctl.clean))


# transpose OTUs
OTU_table.ctl.RA.t <- as.data.frame(t(OTU_table.ctl.RA))

# merge taxonomy and OTU
OTU_table.ctl.RA.tax <- merge(OTU_table.ctl.RA.t, tax_table, by=0)
OTU_table.ctl.RA.tax$Row.names <- NULL

# summarize as top 11 OTU
OTU_table.ctl.RA.sort <- OTU_table.ctl.RA.t[, order(colSums(OTU_table.ctl.RA.t),
                                                    decreasing=TRUE)][1:11,]
OTU_table.ctl.RA.sort.tax <- merge(OTU_table.ctl.RA.sort, tax_table, by=0) # merge by row.name
OTU_table.ctl.RA.sort.tax <- OTU_table.ctl.RA.sort.tax[c("Row.names", "Size", 
                              "kingdom", "phylum", "class", "order", "family", 
                                "genus", "kit_neg_A1", "kit_neg_A2", "kit_neg_B", 
                                  "seq_neg_A1", "seq_neg_A2", "seq_neg_B")]

OTU_table.ctl.RA.sort.tax$Size <- NULL
OTU_table.ctl.RA.sort.tax$kingdom <- NULL

colnames(OTU_table.ctl.RA.sort.tax) <- c("OTU", "Phylum", "Class", "Order", 
                                          "Family", "Genus", "Kit Ai", "kit Aii", 
                                            "kit B", "seq Ai", "seq Aii", "seq B")


# summarize as phyla
OTU_table.ctl.RA.phyla <- aggregate(OTU_table.ctl.RA.tax[, 1:6], 
                                    list(OTU_table.ctl.RA.tax$phylum), sum)

colnames(OTU_table.ctl.RA.phyla) <- c("Phylum", "kit Ai", "kit Aii", "kit B", 
                                        "seq Ai", "seq Aii", "seq B")



# draw tables
pdf("results/tables/Table S1a_NTC_phyla.pdf", height=5, width=20)
grid.table(OTU_table.ctl.RA.phyla, rows=NULL)
dev.off()

pdf("results/tables/Table_S1b_NTC_top10.pdf", height=5, width=20)
grid.table(OTU_table.ctl.RA.sort.tax, rows=NULL)
dev.off()



