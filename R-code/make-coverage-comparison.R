rm(list = ls())


library(plyr)
library(dplyr)
library(ggplot2)

# set your own link to your home repo here, 
# and the rest of the files in the repo should run
homewd = "/Users/carabrook/Developer/mNGS-human-fever"

setwd(homewd)

#load the data for mngs first (the samtools file)

tree 
mngs.reads = read.delim(file = paste0(homewd, "/data/samtools_depth_mngs_hku1.txt"), header = T, stringsAsFactors = F)
head(mngs.reads) # this is a single column with the number of reads per site on the genome
names(mngs.reads) <- "read_depth"
# add column with nucleotide position
mngs.reads$position <- 1: nrow(mngs.reads)

# quick visualization of coverage:

p1 <- ggplot(data=mngs.reads) + geom_line(aes(x=position, y=read_depth))
p1

# nice ! but you really want to plot reads per million, so you need the total number of reads as well:
# load this from the .tsv report

mngs.report <- read.delim(file =paste0(homewd,"/data/report_mngs_hku1.tsv"), header = T, stringsAsFactors = F)
head(mngs.report) #extract total reads and attach to above

mngs.reads$total_reads <- as.numeric(mngs.report$consensus[mngs.report$Assembly=="# total reads"])

#and convert read depth to reads per million
mngs.reads$read_depth_per_million <- mngs.reads$read_depth/(unique(mngs.reads$total_reads)/(1000000))

#and replot as read depth per million

p2 <- ggplot(data=mngs.reads) + geom_line(aes(x=position, y=read_depth_per_million))
p2


#now add identifiers to the data 
mngs.reads$seq_type <- "mNGS"

#and reproduce with the msspe data
msspe.reads = read.delim(file = paste0(homewd, "/data/samtools_depth_msspe_hku1.txt"), header = T, stringsAsFactors = F)
head(msspe.reads) # this is a single column with the number of reads per site on the genome
names(msspe.reads) <- "read_depth"
# add column with nucleotide position
msspe.reads$position <- 1: nrow(msspe.reads)

# quick visualization of coverage:

p3 <- ggplot(data=msspe.reads) + geom_line(aes(x=position, y=read_depth))
p3


msspe.report <- read.delim(file =paste0(homewd,"/data/report_msspe_hku1.tsv"), header = T, stringsAsFactors = F)
head(msspe.report) #extract total reads and attach to above

msspe.reads$total_reads <- as.numeric(msspe.report$consensus[msspe.report$Assembly=="# total reads"])

#and convert read depth to reads per million
msspe.reads$read_depth_per_million <- msspe.reads$read_depth/(unique(msspe.reads$total_reads)/(1000000))


p4 <- ggplot(data=msspe.reads) + geom_line(aes(x=position, y=read_depth_per_million))
p4

#now add identifiers to the data 
msspe.reads$seq_type <- "MSSPE"


#and join
merge.dat <- rbind(mngs.reads, msspe.reads)


p5 <- ggplot(data=merge.dat) + theme_bw() +
      geom_line(aes(x=position, y=read_depth_per_million, color=seq_type))+
      theme(panel.grid = element_blank(), legend.title = element_blank(),
            legend.position = c(.1,.85)) + xlab("genome position") + ylab("read depth per million")
  
p5

#and also compare raw:
p6 <- ggplot(data=merge.dat) + theme_bw() + 
      geom_line(aes(x=position, y=read_depth, color=seq_type)) +
      theme(panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = c(.1,.85)) + xlab("genome position") + ylab("raw read depth")
p6

#and plot together
pout <- cowplot::plot_grid(p6,p5, align = "hv", ncol=1, nrow=2, labels = c("A", "B"), label_size = 22, label_x = c(-.01))

#save to the figures folder
ggsave(file = paste0(homewd,"/figures/mNGS_MSSPE_comparison_read_depth.png"),
       plot = pout,
       units="mm",  
       width=80, 
       height=70, 
       scale=2.5, 
       dpi=300)
