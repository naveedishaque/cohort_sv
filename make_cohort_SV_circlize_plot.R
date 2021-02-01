#

# module load R/3.3.1
# example call for showing all SVs associating with a gene 2+ times, and highlighting those that aossicate 6+ times: 
# R -f make_cohort_SV_circlize_plot.R --args [output from process_SV_for_circos] 6 2 "My title" 

## LOAD LIBS

library(circlize)
library(dplyr)

## PARSE ARGS

min_to_show_sv=1
min_recurrence=5
title="my title";

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
if (length(args)>0) {
  in_file <- args[1]
} 
if (length(args)>1) {
  min_recurrence <- as.numeric(args[2])
} 
if (length(args)>2) {
  min_to_show_sv <- as.numeric(args[3])
} 
if (length(args)>3) {
  title <- args[4]
}

print(paste("file: ", in_file), sep="")
print(paste("min_recurrence: ", min_recurrence), sep="")
print(paste("min_to_show: ", min_to_show_sv), sep="")
title <- paste(title," (minimum shown = ",min_to_show_sv,"; highlight = ",min_recurrence,"+)", sep="")
print(paste(title))

## READ TABLE AND SPLIT

full_table <- read.table(in_file, header=T, sep="\t", fill=T)
head(full_table)
full_table_min <- full_table %>% filter(SCORE >= min_recurrence)
head(full_table_min)

# subset tables

window_table_direct<- full_table %>% filter(TYPE == "windows_SVs_direct")
window_table_near  <- full_table %>% filter(TYPE == "windows_SVs_near")
window_table_close <- full_table %>% filter(TYPE == "windows_SVs_closest")
gene_table_close <- full_table %>% filter(TYPE == "gene_SVs")
gene_table_direct<- full_table %>% filter(TYPE == "gene_SVs_direct")
gene_table_near  <- full_table %>% filter(TYPE == "gene_SVs_near")
del_table        <- full_table %>% filter(TYPE == "del_SVs")
dup_table        <- full_table %>% filter(TYPE == "dup_SVs")

window_table_direct_min  <- full_table_min %>% filter(TYPE == "windows_SVs_direct")
window_table_near_min    <- full_table_min %>% filter(TYPE == "windows_SVs_near")
window_table_close_min   <- full_table_min %>% filter(TYPE == "windows_SVs_closest")

gene_table_close_min <- full_table_min %>% filter(TYPE == "gene_SVs_closest")
gene_table_direct_min<- full_table_min %>% filter(TYPE == "gene_SVs_direct")
gene_table_near_min  <- full_table_min %>% filter(TYPE == "gene_SVs_near")
del_table_min        <- full_table_min %>% filter(TYPE == "del_SVs")
dup_table_min        <- full_table_min %>% filter(TYPE == "dup_SVs")

# create SV bed objects

window_bed_direct_pos1 <- window_table_direct[,c(2,3,3,9)]
window_bed_direct_pos2 <- window_table_direct[,c(5,7,7,9)]
colnames(window_bed_direct_pos1) <- c("chr","start","end","value1")
colnames(window_bed_direct_pos2) <- c("chr","start","end","value1")

window_bed_near_pos1 <- window_table_near[,c(2,3,3,9)]
window_bed_near_pos2 <- window_table_near[,c(5,7,7,9)]
colnames(window_bed_near_pos1) <- c("chr","start","end","value1")
colnames(window_bed_near_pos2) <- c("chr","start","end","value1")

window_bed_close_pos1 <- window_table_close[,c(2,3,3,9)]
window_bed_close_pos2 <- window_table_close[,c(5,7,7,9)]
colnames(window_bed_close_pos1) <- c("chr","start","end","value1")
colnames(window_bed_close_pos2) <- c("chr","start","end","value1")

window_bed_direct_pos1_min <- window_table_direct_min[,c(2,3,3,9)]
window_bed_direct_pos2_min <- window_table_direct_min[,c(5,7,7,9)]
colnames(window_bed_direct_pos1_min) <- c("chr","start","end","value1")
colnames(window_bed_direct_pos2_min) <- c("chr","start","end","value1")

window_bed_near_pos1_min <- window_table_near_min[,c(2,3,3,9)]
window_bed_near_pos2_min <- window_table_near_min[,c(5,7,7,9)]
colnames(window_bed_near_pos1_min) <- c("chr","start","end","value1")
colnames(window_bed_near_pos2_min) <- c("chr","start","end","value1")

window_bed_close_pos1_min <- window_table_close_min[,c(2,3,3,9)]
window_bed_close_pos2_min <- window_table_close_min[,c(5,7,7,9)]
colnames(window_bed_close_pos1_min) <- c("chr","start","end","value1")
colnames(window_bed_close_pos2_min) <- c("chr","start","end","value1")

# reformat gene table

gene_bed_close_min <- gene_table_close_min[,c(2,3,4,11)]
colnames(gene_bed_close_min) <- c("chr","start","end","name")

gene_bed_direct_min <- gene_table_direct_min[,c(2,3,4,11)]
colnames(gene_bed_direct_min) <- c("chr","start","end","name")

gene_bed_near_min <- gene_table_near_min[,c(2,3,4,11)]
colnames(gene_bed_near_min) <- c("chr","start","end","name")

## CHECK TABLES
dim(gene_bed_direct_min)
dim(gene_bed_near_min)
dim(gene_bed_close_min)
dim(window_bed_direct_pos1_min)
dim(window_bed_direct_pos2_min)
dim(window_bed_near_pos1_min)
dim(window_bed_near_pos2_min)
dim(window_bed_close_pos1_min)
dim(window_bed_close_pos2_min)


gene_bed_direct_min <- gene_bed_direct_min[gene_bed_direct_min$start>0,]
gene_bed_near_min <- gene_bed_near_min[gene_bed_near_min$start>0,]
gene_bed_close_min <- gene_bed_close_min[gene_bed_close_min$start>0,]
#window_bed_direct_pos1_min <- window_bed_direct_pos1_min[window_bed_direct_pos1_min$start>0,]
#window_bed_direct_pos2_min <- window_bed_direct_pos2_min[window_bed_direct_pos2_min$start>0,]
#window_bed_near_pos1_min <- window_bed_near_pos1_min[window_bed_near_pos1_min$start>0,]
#window_bed_near_pos2_min <- window_bed_near_pos2_min[window_bed_near_pos2_min$start>0,]
#window_bed_close_pos1_min <- window_bed_close_pos1_min[window_bed_close_pos1_min$start>0,]
#window_bed_close_pos2_min <- window_bed_close_pos2_min[window_bed_close_pos2_min$start>0,]

dim(gene_bed_direct_min)
dim(gene_bed_near_min)
dim(gene_bed_close_min)
nrow(window_bed_direct_pos1_min)
nrow(window_bed_direct_pos2_min)
nrow(window_bed_near_pos1_min)
nrow(window_bed_near_pos2_min)
nrow(window_bed_close_pos1_min)
nrow(window_bed_close_pos2_min)

dim(window_bed_direct_pos1)
dim(window_bed_direct_pos2)
dim(window_bed_near_pos1)
dim(window_bed_near_pos2)
dim(window_bed_close_pos1)
dim(window_bed_close_pos2)

head(window_bed_direct_pos1)
head(window_bed_direct_pos2)
tail(window_bed_direct_pos1)
tail(window_bed_direct_pos2)

## PRINT

outfile = paste(args[1],".SVcohort.low",min_to_show_sv,".high",min_recurrence,".pdf", sep="")

pdf(outfile, width=8, height=8)

circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
title(main=paste(title," (SV direct)", sep=""))
circos.genomicLabels(gene_bed_direct_min, labels.column = 4, side = "outside",col="black")
circos.genomicIdeogram()
circos.track(ylim=c(0, 1), panel.fun= function(x,y) {circos.text(CELL_META$xcenter, CELL_META$ycenter,gsub("chr","",CELL_META$sector.index))}, track.height=0.05, bg.border=NA)
circos.genomicLink(window_bed_direct_pos1    , window_bed_direct_pos1    , col="#FFFFFF", border=0.5, lwd=2)
circos.genomicLink(window_bed_direct_pos1_min, window_bed_direct_pos2_min, col="#FF000030", border=0.5, lwd=2)

circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
title(main=paste(title," (SV near)", sep=""))
circos.genomicLabels(gene_bed_near_min, labels.column = 4, side = "outside",col="black")
circos.genomicIdeogram()
circos.track(ylim=c(0, 1), panel.fun= function(x,y) {circos.text(CELL_META$xcenter, CELL_META$ycenter,gsub("chr","",CELL_META$sector.index))}, track.height=0.05, bg.border=NA)
circos.genomicLink(window_bed_near_pos1    , window_bed_near_pos1    , col="#FFFFFF", border=0.5, lwd=2)
circos.genomicLink(window_bed_near_pos1_min, window_bed_near_pos2_min, col="#FF000030", border=0.5, lwd=2)

circos.clear()
circos.initializeWithIdeogram(plotType = NULL)
title(main=paste(title," (SV closest)", sep=""))
circos.genomicLabels(gene_bed_close_min, labels.column = 4, side = "outside",col="black")
circos.genomicIdeogram()
circos.track(ylim=c(0, 1), panel.fun= function(x,y) {circos.text(CELL_META$xcenter, CELL_META$ycenter,gsub("chr","",CELL_META$sector.index))}, track.height=0.05, bg.border=NA)
circos.genomicLink(window_bed_close_pos1    , window_bed_close_pos1    , col="#FFFFFF", border=0.5, lwd=2)
circos.genomicLink(window_bed_close_pos1_min, window_bed_close_pos2_min, col="#FF000030", border=0.5, lwd=2)

dev.off()
