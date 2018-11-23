##############################################################################
#
# Figure 3
# OTU-level stack barcharts
#
# Adapted from B Torondel, JHJ Ensink, O Gundogdu, UZ Ijaz, 
# J Parkhill, F Abdelahi, V-A Nguyen, S Sudgen, W Gibson, 
# AW Walker, and C Quince.
# Assessment of the influence of intrinsic environmental and 
# geographical factors on #the bacterial ecology of pit latrines
# Microbial Biotechnology, 9(2):209-223, 2016. DOI:10.1111/1751-7915.12334 
#
# Code kindly made available here: http://userweb.eng.gla.ac.uk/umer.ijaz)
#
##############################################################################

library("readxl")
library("tidyr")
library("ggplot2")
library("dplyr")
library("tibble")

# read in data
OTU_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared", header = TRUE, row.names = 2)
OTU_table <- OTU_table[-c(1,2)] #tidy

shallow.community.remove <- c("il_cmp_d2_2", "il_cmpphg_d4_3")
OTU_table <- OTU_table[!(row.names(OTU_table) %in% shallow.community.remove), ]

tax_table <- read.table("data/mothur/phage_ecol.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy", header=TRUE)
tax_table <- separate(tax_table, Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), ";")

genus.ID <- t(data.frame(row.names = tax_table$OTU, paste(tax_table$OTU, sep = "; ", tax_table$genus)))
OTU_table.genus <- rbind(OTU_table, genus.ID)

colnames(OTU_table.genus) <- OTU_table.genus[103, ]
OTU_table.genus <- OTU_table.genus[-103, ]

OTU_table.genus <- data.frame(lapply(OTU_table.genus, function(x) as.numeric(as.character(x))),check.names=F, row.names = rownames(OTU_table.genus))


# get grouping information
grouping_info <- data.frame(row.names=row.names(OTU_table), row.names(OTU_table))
grouping_info <- separate(grouping_info, row.names.OTU_table., c("site", "treatment", "time", "replicate"), "_")
grouping_info$day <- ifelse(grouping_info$time=="d1", "1 dpt", ifelse(grouping_info$time=="d2", "2 dpt", ifelse(grouping_info$time=="d3", "3 dpt", ifelse(grouping_info$time=="d4", "4 dpt", ifelse(grouping_info$time=="d5", "5 dpt", NA)))))

grouping_info.ca <- grouping_info[grouping_info$site=="ca" & grouping_info$treatment!="ctl",]
grouping_info.il <- grouping_info[grouping_info$site=="il" & grouping_info$treatment!="ctl",]

## ceca ##
OTU_table.ca <- rownames_to_column(OTU_table.genus)
OTU_table.ca <- OTU_table.ca  %>% filter(rowname %in% row.names(grouping_info.ca)) 
OTU_table.ca <- column_to_rownames(OTU_table.ca, "rowname")


# apply proportion normalisation
x.ca <- 100*(OTU_table.ca/rowSums(OTU_table.ca))
x.ca <- x.ca[, order(colSums(x.ca), decreasing=TRUE)]

# extract list of top N Taxa
N.ca <- 11
taxa_list.ca <- colnames(x.ca)[1:N.ca]
N.ca <- length(taxa_list.ca)

# generate a new table with everything added to Others
new_x.ca <- data.frame(x.ca[,colnames(x.ca) %in% taxa_list.ca], Others=rowSums(x.ca[,!colnames(x.ca) %in% taxa_list.ca]))


# change the Type=grouping_info[,1] should you desire any other grouping of panels
df.ca <- NULL
for (i in 1:dim(new_x.ca)[2]){
  tmp <- data.frame(row.names=NULL, Sample=rownames(new_x.ca), 
                    Taxa = rep(colnames(new_x.ca)[i], dim(new_x.ca)[1]), Value=new_x.ca[,i], 
                    Type= grouping_info.ca[,5])
  if(i==1){df.ca <- tmp} else {df.ca <- rbind(df.ca,tmp)}
}

df.ca$taxa_format <- gsub(".100.", " (100)", df.ca$Taxa)
df.ca$taxa_format <- gsub(".98.", " (98)", df.ca$taxa_format)



colours.ca <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#C20088","#9DCC00","#003380")


p.ca <- ggplot(df.ca, aes(Sample, Value, fill=taxa_format)) + 
          geom_bar(stat="identity") + facet_grid(. ~ Type, drop=TRUE, scale="free", space="free_x") +
        scale_fill_manual(values=colours.ca[1:(N.ca+1)], 
          breaks=c("Otu0001..Lactobacillus (100)", "Otu0003..Enterobacteriaceae_unclassified (100)", 
            "Otu0004..Faecalibacterium (100)", "Otu0006..Clostridium_IV (98)", 
            "Otu0008..Ruminococcaceae_unclassified (100)", "Otu0009..Clostridiales_unclassified (100)", 
            "Otu0010..Lachnospiraceae_unclassified (100)", "Otu0012..Lachnospiraceae_unclassified (100)", 
            "Otu0013..Campylobacter (100)", "Otu0014..Lachnospiraceae_unclassified (100)", 
            "Otu0015..Lachnospiraceae_unclassified (100)", "Others")) +
          ylab("Relative abundance (%)") +
          labs(fill = "Genus") +
          scale_y_continuous(expand = c(0,0)) + 
          scale_x_discrete(labels=c("ca_cmp_d1_1" = "Cj_1", "ca_cmp_d1_2" = "Cj_2", "ca_cmp_d1_3" = "Cj_3", 
              "ca_cmp_d1_4" = "Cj_4", "ca_cmp_d1_5" = "Cj_5", "ca_cmp_d2_1" = "Cj_1", "ca_cmp_d2_2" = "Cj_2",
              "ca_cmp_d2_3" = "Cj_3", "ca_cmp_d2_4" = "Cj_4", "ca_cmp_d2_5" = "Cj_5", "ca_cmp_d3_1" = "Cj_1", 
              "ca_cmp_d3_2" = "Cj_2", "ca_cmp_d3_3" = "Cj_3", "ca_cmp_d3_4" = "Cj_4", "ca_cmp_d3_5" = "Cj_5", 
              "ca_cmp_d4_1" = "Cj_1", "ca_cmp_d4_2" = "Cj_2", "ca_cmp_d4_3" = "Cj_3", "ca_cmp_d4_4" = "Cj_4", 
              "ca_cmp_d4_5" = "Cj_5", "ca_cmp_d5_1" = "Cj_1", "ca_cmp_d5_2" = "Cj_2", "ca_cmp_d5_3" = "Cj_3", 
              "ca_cmp_d5_4" = "Cj_4", "ca_cmpphg_d1_1" = "Cj_phg_1", "ca_cmpphg_d1_2" = "Cj_phg_2", 
              "ca_cmpphg_d1_3" = "Cj_phg_3", "ca_cmpphg_d1_4" = "Cj_phg_4", "ca_cmpphg_d1_5" = "Cj_phg_5", 
              "ca_cmpphg_d2_1" = "Cj_phg_1", "ca_cmpphg_d2_2" = "Cj_phg_2", "ca_cmpphg_d2_3" = "Cj_phg_3", 
              "ca_cmpphg_d2_4" = "Cj_phg_4", "ca_cmpphg_d2_5" = "Cj_phg_5", "ca_cmpphg_d3_1" = "Cj_phg_1", 
              "ca_cmpphg_d3_2" = "Cj_phg_2", "ca_cmpphg_d3_3" = "Cj_phg_3", "ca_cmpphg_d3_4" = "Cj_phg_4", 
              "ca_cmpphg_d3_5" = "Cj_phg_5", "ca_cmpphg_d4_1" = "Cj_phg_1", "ca_cmpphg_d4_2" = "Cj_phg_2", 
              "ca_cmpphg_d4_3" = "Cj_phg_3", "ca_cmpphg_d4_4" = "Cj_phg_4", "ca_cmpphg_d4_5" = "Cj_phg_5", 
              "ca_cmpphg_d5_1" = "Cj_phg_1", "ca_cmpphg_d5_2" = "Cj_phg_2", "ca_cmpphg_d5_3" = "Cj_phg_3", 
              "ca_cmpphg_d5_4" = "Cj_phg_4", "ca_cmpphg_d5_5" = "Cj_phg_5")) +
          theme_bw() + 
          theme(strip.background = element_rect(fill="white"), panel.spacing = unit(0.3, "lines"), 
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8, colour = "black"), 
                axis.text.y=element_text(colour = "black"))


tiff(filename="results/figures/fig3b_ceca.tiff", width=250, height=100, units="mm", res=300)
p.ca
dev.off()



## ileum ##
OTU_table.il <- rownames_to_column(OTU_table.genus)
OTU_table.il <- OTU_table.il  %>% filter(rowname %in% row.names(grouping_info.il)) 
OTU_table.il <- column_to_rownames(OTU_table.il, "rowname")

# apply proportion normalisation
x.il <- 100*(OTU_table.il/rowSums(OTU_table.il))
x.il <- x.il[, order(colSums(x.il), decreasing=TRUE)]

# extract list of top N Taxa
N.il <- 11
taxa_list.il <- colnames(x.il)[1:N.il]
N.il <- length(taxa_list.il)

# generate a new table with everything added to Others
new_x.il <- data.frame(x.il[,colnames(x.il) %in% taxa_list.il], Others=rowSums(x.il[,!colnames(x.il) %in% taxa_list.il]))


# change the Type=grouping_info[,1] should you desire any other grouping of panels
df.il <- NULL
for (i in 1:dim(new_x.il)[2]){
  tmp <- data.frame(row.names=NULL, Sample=rownames(new_x.il), 
                    Taxa = rep(colnames(new_x.il)[i], dim(new_x.il)[1]), Value=new_x.il[,i], 
                    Type= grouping_info.il[,5])
  if(i==1){df.il <- tmp} else {df.il <- rbind(df.il, tmp)}
}

df.il$taxa_format <- gsub(".100.", " (100)", df.il$Taxa)
df.il$taxa_format <- gsub(".98.", " (98)", df.il$taxa_format)
df.il$taxa_format <- gsub(".93.", " (93)", df.il$taxa_format)
df.il$taxa_format <- gsub(".91.", " (91)", df.il$taxa_format)

# reargange so that they have the same scheme as caecal data
#colours.il <- c("#F0A3FF", "#993F00", "#2BCE48",  "#808080", "#8F7C00", "#C20088", "#FFFF00", "#FFA8BB", "#FF0010", "#00998F", "#740AFF", "#990000")
colours.il <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#C20088", "#8F7C00","#9DCC00", "#94FFB5", "#003380")



p.il <- ggplot(df.il, aes(Sample, Value, fill=taxa_format)) + 
          geom_bar(stat="identity") + facet_grid(. ~ Type, drop=TRUE, scale="free", space="free_x") +
          scale_fill_manual(values=colours.il[1:(N.il+1)], 
            breaks=c("Otu0001..Lactobacillus (100)", "Otu0002..Lactobacillus (100)", 
              "Otu0003..Enterobacteriaceae_unclassified (100)", "Otu0005..Romboutsia (100)", 
              "Otu0007..Clostridiales_unclassified (100)", "Otu0011..Prevotella (100)", 
              "Otu0013..Campylobacter (100)", "Otu0016..Enterococcus (100)", 
              "Otu0044..Enterobacteriaceae_unclassified (91)", "Otu0040..Lactobacillus (100)", 
              "Otu0064..Clostridium_sensu_stricto (100)", "Others")) +
          ylab("Relative abundance (%)") +
          labs(fill = "Genus") +
          scale_y_continuous(expand = c(0,0)) + 
          scale_x_discrete(labels=c("il_cmp_d2_1" = "Cj_1", "il_cmp_d2_2" = "Cj_2", "il_cmp_d2_3" = "Cj_3", 
              "il_cmp_d2_4" = "Cj_4", "il_cmp_d2_5" = "Cj_5", "il_cmp_d3_1" = "Cj_1", "il_cmp_d3_2" = "Cj_2",
              "il_cmp_d3_3" = "Cj_3", "il_cmp_d3_4" = "Cj_4", "il_cmp_d3_5" = "Cj_5", "il_cmp_d4_1" = "Cj_1",
              "il_cmp_d4_2" = "Cj_2", "il_cmp_d4_3" = "Cj_3", "il_cmp_d4_4" = "Cj_4", "il_cmp_d4_5" = "Cj_5",
              "il_cmp_d5_1" = "Cj_1", "il_cmp_d5_2" = "Cj_2", "il_cmp_d5_3" = "Cj_3", "il_cmp_d5_4" = "Cj_4",
              "il_cmp_d5_5" = "Cj_5", "il_cmpphg_d2_1" = "Cj_phg_1", "il_cmpphg_d2_2" = "CJ_phg_2", 
              "il_cmpphg_d2_3" = "Cj_phg_3", "il_cmpphg_d2_4" = "Cj_phg_4", "il_cmpphg_d2_5" = "Cj_phg_5",
              "il_cmpphg_d3_1" = "Cj_phg_1", "il_cmpphg_d3_2" = "Cj_phg_2", "il_cmpphg_d3_3" = "Cj_phg_3",
              "il_cmpphg_d3_4" = "Cj_phg_4", "il_cmpphg_d3_5" = "Cj_phg_5", "il_cmpphg_d4_1" = "Cj_phg_1",
              "il_cmpphg_d4_2" = "Cj_phg_2", "il_cmpphg_d4_3" = "Cj_phg_3", "il_cmpphg_d4_4" = "Cj_phg_4",
              "il_cmpphg_d4_5" = "Cj_phg_5", "il_cmpphg_d5_1" = "Cj_phg_1", "il_cmpphg_d5_2" = "Cj_phg_2",
              "il_cmpphg_d5_3" = "Cj_phg_3", "il_cmpphg_d5_4" = "Cj_phg_4", "il_cmpphg_d5_5" = "Cj_phg_5")) +
          theme_bw() + 
          theme(strip.background = element_rect(fill="white"), panel.spacing = unit(0.3, "lines"), 
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8, colour = "black"), 
          axis.text.y=element_text(colour = "black"))

tiff(filename="results/figures/fig3a_ileum.tiff", width=250, height=100, units="mm", res=300)
p.il
dev.off()

