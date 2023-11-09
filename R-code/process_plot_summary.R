rm(list = ls())

library(plyr)
library(dplyr)
library(ggplot2)

# set your own link to your home repo here, 
# and the rest of the files in the repo should run
homewd = "/Users/carabrook/Developer/mNGS-human-fever"

setwd(homewd)

#load the data
sum.dat = read.csv(file = paste0(homewd,"/data/gce_sample_summary.csv"), header=T, stringsAsFactors = F)
head(sum.dat)
names(sum.dat)

sum.dat = arrange(sum.dat, total_reads)
sum.dat$sample_number = seq(1, length(sum.dat$sample_name),1)
sum.dat$sample_number = factor(sum.dat$sample_number, levels=unique(sum.dat$sample_number))

#visualize total rea
p1 <- ggplot(data=sum.dat) + geom_bar(aes(x=sample_number, y = log(total_reads), fill=sample_type), stat="identity") + theme_bw() + 
  theme(legend.title = element_blank(), panel.grid = element_blank(),  axis.title.x = element_blank(), legend.position=c(.8,.21), legend.text = element_text(size=8), axis.text.y = element_text(size=16), axis.title.y = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("total reads") + 
  scale_y_discrete(limits=c(2.302585, 4.60517, 6.907755, 9.21034, 11.51293, 13.81551, 16.1181), labels=c( "10", "100", "1000",  "10000", "100000", "1000000", "10000000"))
print(p1)

#save to the figures folder
ggsave(file = paste0(homewd,"/figures/GCE_QC_by_sample_type.png"),
       plot = p1,
       units="mm",  
       width=80, 
       height=50, 
       scale=3, 
       dpi=300)



# How many reads are remaining after filter by CZID?
# You can see how only the NP swab samples were of high
# quality to be retained
sum.dat = arrange(sum.dat, reads_after_cdhitdup)
sum.dat$sample_number = seq(1, length(sum.dat$sample_name),1)
sum.dat$sample_number = factor(sum.dat$sample_number, levels=unique(sum.dat$sample_number))

p2 <- ggplot(data=sum.dat) +
  geom_bar(aes(x=sample_number, y = log(reads_after_cdhitdup), fill=sample_type), stat="identity", show.legend = F) +
  theme_bw() + theme(axis.text.y = element_text(size=16), axis.title.y = element_text(size=18), panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("reads after filtering") +
  scale_y_discrete(limits=c(2.302585, 4.60517, 6.907755, 9.21034, 11.51293, 13.81551, 16.1181), labels=c( "10", "100", "1000",  "10000", "100000", "1000000", "10000000"))
print(p2)


p.out1 <- cowplot::plot_grid(p1,p2, nrow=1, ncol=2)

print(p.out1)

ggsave(file = paste0(homewd,"/figures/QC_GCE_by_sample_type.png"),
       plot=p.out1,
       units="mm",  
       width=100, 
       height=40, 
       scale=3, 
       dpi=300)


# Now get proportion of total reads remaining after each filtering step in CZID
# Note - I wrote this script many years ago, so it could 
# definitely be formatted more efficiently -- but it works

names(sum.dat)

# Here, just selecting the columns of interest -- the different numbers of reads remaining after
# each of the curation/filtration steps that are run by czid. So you can see where you lose reads
# and what proportion of the original are remaining at the end of the curation.
prop.reads <- dplyr::select(sum.dat, sample_name, sample_type, total_reads, nonhost_reads, reads_after_star, reads_after_trimmomatic, reads_after_priceseq, reads_after_cdhitdup)
head(prop.reads)
names(prop.reads) <- c("sample_name", "sample_type", "total_reads", "nonhost-reads", "filter-star", "filter-trimmomatic", "filter-priceseq", "filter-cdhitup")
head(prop.reads)

# First, pull just the total reads
total.reads <- dplyr::select(prop.reads, sample_name, sample_type, total_reads)
head(total.reads)

# Drop from other file
prop.reads <- dplyr::select(prop.reads, -(total_reads))

# Now, melt this into long format by sample name and type
library(reshape2)
reads.long <- melt(prop.reads, id.vars = c("sample_name", "sample_type"), variable.name = "read_type", value.name = "reads_remaining")
head(reads.long)

#and merge with total reads
reads.long <- merge(reads.long, total.reads, by=c("sample_name", "sample_type"), all.x = T)
head(reads.long)

# and get proportion of the total original that are lost with each step
# order of czid is star - trimmomatic - priceseq - cdhit - then several small steps (see poster)

# split by sample and run function
sample.split <- dlply(reads.long,.(sample_name))

read.lost <- function(df){
  
  
  df$reads_lost <- NA
  df$reads_lost[df$read_type=="filter-star"] <- df$total_reads[df$read_type=="filter-star"] - df$reads_remaining[df$read_type=="filter-star"] 
  df$reads_lost[df$read_type=="filter-trimmomatic"] <- df$reads_remaining[df$read_type=="filter-star"] - df$reads_remaining[df$read_type=="filter-trimmomatic"] 
  df$reads_lost[df$read_type=="filter-priceseq"] <-df$reads_remaining[df$read_type=="filter-trimmomatic"]  - df$reads_remaining[df$read_type=="filter-priceseq"] 
  df$reads_lost[df$read_type=="filter-cdhitup"] <-df$reads_remaining[df$read_type=="filter-priceseq"]  - df$reads_remaining[df$read_type=="filter-cdhitup"] 
  
  
  #then, the remaining difference between the current total and the non-host reads are other filtering steps
  
  df$reads_lost[df$read_type=="nonhost-reads"] <-df$reads_remaining[df$read_type=="filter-cdhitup"]  - df$reads_remaining[df$read_type=="nonhost-reads"] 
  
  #add in the total retained
  df.1 <- df[df$read_type=="nonhost-reads",]
  df.1$reads_lost=df.1$reads_remaining

  df$read_type <- as.character(df$read_type)
  df$read_type[df$read_type=="nonhost-reads"] <- "filter-final"
  
  df <- rbind(df, df.1)
  #and calculate the proportion of total cut out at each step - and proportion remaining
  df$proportion_lost <- df$reads_lost/df$total_reads # these should all add to 1 if summed
  
  #for sorting, add proportion retained
  df$prop_retained <- df$reads_remaining[df$read_type=="nonhost-reads"]/unique(df$total_reads)
  
  return(df)

}

#and run function
reads.lost.df <-data.table::rbindlist(lapply(sample.split, read.lost))
head(reads.lost.df) # proportion lost at each step should total to 1 for each sample

# and gut.check -- do the proportions add to 1?
check.tot <- ddply(reads.lost.df, .(sample_name), summarise, total_prop = sum(proportion_lost))
check.tot # pretty good!
# and plot - these should add

reads.lost.df$read_type <- factor(reads.lost.df$read_type, levels=c("filter-star", "filter-trimmomatic", "filter-priceseq", "filter-cdhitup", "filter-final", "nonhost-reads"))
reads.lost.df$sample_type <- factor(reads.lost.df$sample_type, levels=c("water", "whole blood", "plasma", "serum", "nasopharyngeal swab"))
#arrange by those with the highest retained
reads.lost.df <- arrange(reads.lost.df, sample_type,  prop_retained, read_type)
reads.lost.df$sample_name <- factor(reads.lost.df$sample_name, levels = c(unique(reads.lost.df$sample_name)))
tmp = scales::hue_pal()(length(c(unique(reads.lost.df$read_type), unique(reads.lost.df$read_type) )))
colz = tmp[1:length(unique(reads.lost.df$read_type))]
fillz = tmp[(length(unique(reads.lost.df$read_type))+1):length(tmp)]

p.out2 <- ggplot(data=reads.lost.df) + 
  geom_bar(aes(x=sample_name, y=proportion_lost, 
               fill=read_type, color=sample_type), stat="identity") +
  scale_color_manual(values=colz) + scale_fill_manual(values=fillz) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=18), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        legend.position = "top", legend.title=element_blank()) + ylab("proportion of reads")
print(p.out2)


ggsave(file = paste0(homewd,"/figures/QC_GCE_prop_reads_spacer.png"),
       plot = p.out2,
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)

p.out3 <- ggplot(data=reads.lost.df) + 
  geom_bar(aes(x=sample_name, y=proportion_lost, 
               fill=read_type), stat="identity") +
  scale_fill_manual(values=fillz) +#scale_color_manual(values=colz) + 
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=18), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        legend.position = "top", legend.title=element_blank()) + ylab("proportion of reads")
print(p.out3)

ggsave(file = paste0(homewd,"/figures/QC_GCE_prop_reads.png"),
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)


#and look at some QC metrics by sample type: QC, compression ratio, non-host-reads
head(sum.dat)
sum.dat$sample_type = factor(sum.dat$sample_type, levels = c("water", "whole blood", "plasma", "serum", "nasopharyngeal swab"))

p4 <- ggplot(sum.dat) + geom_violin(aes(x=sample_type, y=quality_control, fill=sample_type), draw_quantiles=c(0.25, 0.5, 0.75), show.legend = F) +
  coord_flip() + theme_bw() + theme(panel.grid = element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=18), axis.text = element_text(size=14)) + ylab("quality control") 
print(p4)

p5 <- ggplot(sum.dat) + geom_violin(aes(x=sample_type, y=compression_ratio, fill=sample_type), draw_quantiles=c(0.25, 0.5, 0.75), show.legend = F) +
  coord_flip() + theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=18), axis.text = element_text(size=14)) + ylab("compression ratio") 
print(p5)

p6 <- ggplot(sum.dat) + geom_violin(aes(x=sample_type, y=log(nonhost_reads), fill=sample_type), draw_quantiles=c(0.25, 0.5, 0.75), show.legend = F) +
  coord_flip() + theme_bw() + theme(panel.grid = element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=18), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text = element_text(size=14)) + ylab("log(non-host reads)") 
print(p6)

p.out4 <- cowplot::plot_grid(p4,p5,p6, nrow = 1, ncol=3, rel_widths = c(1.5,1,1))
print(p.out4)
ggsave(file = paste0(homewd,"/figures/QC_GCE_by_sample_type.png"),
       plot = p.out4,
       units="mm",  
       width=100, 
       height=50, 
       scale=3, 
       dpi=300)

head(sum.dat)

#just a gut-check: we spike in the ERCC to track the sequencing process - 
#we should see these positively correlated such that we are not losing reads
#somewhere in the pipeline.
ggplot(data=sum.dat) + geom_point(aes(x=log(total_reads), y=log(total_ercc_reads), color=sample_type))

# Next, get proportion of sample of each type of microbe (GCE human_analyses script)

