rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(beeswarm)


# set your own link to your home repo here, 
# and the rest of the files in the repo should run
homewd = "/Users/carabrook/Developer/mNGS-human-fever"

setwd(homewd)


setwd("/Users/caraebrook/Documents/R/R_repositories/CZB/GCE_heatmap/")
#load the output data from the human GCE - and build into a master dataset
heat.NP.virus <- read.csv(file="heatmap_NPswabs_virus.csv", header = T, stringsAsFactors = F)
heat.NP.eukaryota <- read.csv(file="heatmap_NPswabs_eukaryota.csv", header = T, stringsAsFactors = F)
heat.NP.bacteria <- read.csv(file="heatmap_NPswabs_bacteria.csv", header = T, stringsAsFactors = F)
heat.NP.virus$taxon_type <- "virus"
heat.NP.eukaryota$taxon_type <- "eukaryota"
heat.NP.bacteria$taxon_type <- "bacteria"
heat.NP <- rbind(heat.NP.virus, heat.NP.bacteria, heat.NP.eukaryota)

heat.SP.virus <- read.csv(file="heatmap_serum_plasma_virus.csv", header = T, stringsAsFactors = F)
heat.SP.eukaryota <- read.csv(file="heatmap_serum_plasma_eukaryota.csv", header = T, stringsAsFactors = F)
heat.SP.bacteria <- read.csv(file="heatmap_serum_plasma_bacteria.csv", header = T, stringsAsFactors = F)
heat.SP.virus$taxon_type <- "virus"
heat.SP.eukaryota$taxon_type <- "eukaryota"
heat.SP.bacteria$taxon_type <- "bacteria"
heat.SP <- rbind(heat.SP.virus, heat.SP.bacteria, heat.SP.eukaryota)


heat.WB.virus <- read.csv(file="heatmap_wholeblood_virus.csv", header = T, stringsAsFactors = F)
heat.WB.eukaryota <- read.csv(file="heatmap_wholeblood_eukaryota.csv", header = T, stringsAsFactors = F)
heat.WB.bacteria <- read.csv(file="heatmap_wholeblood_bacteria.csv", header = T, stringsAsFactors = F)
heat.WB.virus$taxon_type <- "virus"
heat.WB.eukaryota$taxon_type <- "eukaryota"
heat.WB.bacteria$taxon_type <- "bacteria"

heat.WB <- rbind(heat.WB.virus, heat.WB.bacteria, heat.WB.eukaryota)


#and combine
heat.dat <- rbind(heat.NP, heat.SP, heat.WB)
#heat.dat <- read.csv(file="IDseq_GCE_heatmap.csv", header = T, stringsAsFactors = F)
head(heat.dat)
length(unique(heat.dat$sample_name)) #313. now 309 must be ignoring the weird pathogen groups


#now pair it with most up-to-date appropriate metadata
meta.dat <- read.csv(file="/Users/caraebrook/Documents/R/R_repositories/CZB/NextSeq_dat_for_upload.csv", header=T, stringsAsFactors = F )
#dat <- read.csv(file = "GCE_plate1_for_upload.csv", header = T, stringsAsFactors = F)
#head(dat)
head(meta.dat)
length(unique(meta.dat$Sample.Name)) #323...what are the 7 that did not report sequences?
names(meta.dat)[1] <- "sample_name"

#setdiff(unique(meta.dat$sample_name), unique(heat.dat$sample_name))
#"RR034H_157_wbldRa_S155", "RR034H_163_wbldRa_S161", "RR034H_164_wbldRa_S162", "RR034H_188_wbldRa_S186", "RR034H_199_wbldRa_S197", "RR034H_201_wbldRa_S199", "RR034H_304_plRa_S302" 
#what are the corresponting host IDs? (from before)
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_157_wbldRa_S155"] #004-FAR
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_163_wbldRa_S161"] #006-TSL
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_164_wbldRa_S162"] #011-MJR
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_188_wbldRa_S186"] #002-BHK
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_199_wbldRa_S197"] #002-TLR
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_201_wbldRa_S199"] #004-TLR
#meta.dat$Host.ID[meta.dat$sample_name=="RR034H_304_plRa_S302"] #03091-16

#these are scattered all over...maybe just failed libraries?


#and merge
merge.dat <- merge(heat.dat, meta.dat,  by="sample_name")
head(merge.dat)
length(unique(merge.dat$sample_name)) #313.

#cull the duplicate waters
unique(merge.dat[duplicated(merge.dat),]$Water.Control) #all duplicates are water
merge.dat <- merge.dat[!duplicated(merge.dat),]

#now do some analysis by metadata...?
#can you reconstruct the heat map?
#total reads >5 and would like to include bp length >50 only...
heat.sub <- subset(merge.dat, NT_r>10)
head(heat.sub)
#for  viruses, you want to plot total reads. for eukaryota, plot the reads per million
heat.sub$Sample.Type = factor(heat.sub$Sample.Type, levels = c("nasopharyngeal swab", "serum", "whole blood", "plasma", "water"))
vir.sub = subset(heat.sub, taxon_type=="virus")
head(vir.sub)
p1 <- ggplot(data=vir.sub) + geom_tile(aes(x=sample_name, y= taxon_name, fill=log(NT_r))) + 
  scale_fill_gradient(low = "white", high = "red") + facet_grid(Sample.Type~.)
print(p1)

#and just NP
NP.sub <- subset(heat.sub, Sample.Type=="nasopharyngeal swab" | Sample.Type=="water")
NP.sub$taxon_type = factor(NP.sub$taxon_type, levels = c("bacteria", "eukaryota", "virus"))
NP.sub <- arrange(NP.sub, taxon_type)
NP.sub$taxon_name = factor(NP.sub$taxon_name, levels = unique(NP.sub$taxon_name))
unique(NP.sub$Primary.Diagnosis)
unique(NP.sub$Admission.Type)
NP.sub$Admission.Type[is.na(NP.sub$Admission.Type)] <- "not admitted"
NP.sub$Admission.Type = factor(NP.sub$Admission.Type, levels = c("not admitted", "fever+cough"))
NP.sub <- arrange(NP.sub,  Admission.Type)
head(NP.sub)
NP.sub$sampleID <- NA
for (i in 1:length(unique(NP.sub$sample_name[NP.sub$Admission.Type=="not admitted" & NP.sub$Water.Control==" No"]))){
  NP.sub$sampleID[NP.sub$sample_name==unique(NP.sub$sample_name[NP.sub$Admission.Type=="not admitted" & NP.sub$Water.Control==" No"])[i]] <- paste0("control-healthy-", i)
}

for (i in 1:length(unique(NP.sub$sample_name[NP.sub$Admission.Type=="not admitted" & NP.sub$Water.Control=="Yes"]))){
  NP.sub$sampleID[NP.sub$sample_name==unique(NP.sub$sample_name[NP.sub$Admission.Type=="not admitted" & NP.sub$Water.Control=="Yes"])[i]] <- paste0("control-water-", i)
}


for (i in 1:length(unique(NP.sub$sample_name[NP.sub$Admission.Type=="fever+cough"]))){
  NP.sub$sampleID[NP.sub$sample_name==unique(NP.sub$sample_name[NP.sub$Admission.Type=="fever+cough"])[i]] <- paste0("fever-cough-", i)
}

#NP.sub <- arrange(NP.sub, sampleID)
NP.sub$sampleID = factor(NP.sub$sampleID, levels=unique(NP.sub$sampleID))


#virus - protein coding reads
p2 <- ggplot(data=subset(NP.sub, taxon_type=="virus")) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NR_rpm))) + 
  scale_fill_gradient2(low = "white", mid="yellow", high = "red") 
print(p2)

sub.dat = subset(NP.sub, taxon_type=="virus" & Sample.Type=="nasopharyngeal swab")
sum.tot = ddply(sub.dat, .(taxon_name), summarise, tot=length(taxon_name))
sum.tot = arrange(sum.tot, desc(tot))
sum.tot$taxon_name = factor(sum.tot$taxon_name, levels= unique(sum.tot$taxon_name))


sub.dat$taxon_name = factor(sub.dat$taxon_name, levels = unique(sum.tot$taxon_name))


sub.dat = arrange(sub.dat, taxon_name)
p2 <- ggplot(data=sub.dat) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NR_rpm))) + 
  scale_fill_gradient2(low = "white", mid="yellow", high = "red") +theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), # axis.text.y = element_text(size=16), 
        axis.text= element_blank(), 
        legend.text = element_text(size=16),
        axis.ticks.y=element_blank(), legend.position = c(.15,.85),
        plot.margin = unit(c(.5,.5,2.9,0), "lines"))
print(p2)

#and look at prevalence

sum.tot$N <- sum(sum.tot$tot)
sum.tot$prev <- sum.tot$tot/sum.tot$N

p3 <- ggplot(data=sum.tot) + geom_bar(aes(x=taxon_name, y=prev, fill=prev), stat = "identity") + 
  theme_bw() + coord_flip() + #scale_y_discrete(position = "top") +
  scale_fill_gradient2(low = "white", mid="yellow", high = "red") +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size=22), 
        axis.text = element_text( size=22), 
        plot.margin = unit(c(.5,0,.5,.5), "lines"),
        legend.text = element_text(size=16),
        legend.position = c(.85,.85), legend.title = element_blank()) + ylab("prevalence")
print(p3)

pout <- cowplot::plot_grid(p3, p2, nrow=1, ncol=2, rel_widths = c(1.5,1) )
pout

ggsave(file = "GCE_NP_swab_viruses.pdf",
       units="mm",  
       width=140, 
       height=80, 
       scale=3, 
       dpi=300)




#compare all - non-protein mapping
p3 <- ggplot(data=NP.sub) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + 
  scale_fill_gradient2(low = "white", mid="yellow", high = "red")  + facet_grid(~taxon_type)
print(p3)

#but not as facets
p4 <- ggplot(data=NP.sub) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + 
  scale_fill_gradient2(low = "yellow", mid="orange", high = "red") # + facet_grid(~taxon_type)
print(p4)

#then flip everything
p5 <- ggplot(data=NP.sub) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + theme_bw() + 
  theme(axis.title = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=90, size=4.5), panel.grid = element_blank(), 
        axis.text.y=element_text(size=3.3), strip.background = element_blank(), panel.spacing=unit(0, "lines"), legend.text = element_text(size=5)) +
  scale_fill_gradient(low ="yellow",  high = "red") + scale_x_discrete(position = "top") #+ facet_grid(taxon_type~.)
print(p5)

#and save
ggsave(file = "GCE_NP_swab_all_pathogens.pdf",
       units="mm",  
       width=70, 
       height=150, 
       scale=3, 
       dpi=300)


#and do for the other sample types
WB.sub <- subset(heat.sub, Sample.Type=="whole blood" | Sample.Type=="water")
WB.sub$taxon_type = factor(WB.sub$taxon_type, levels = c("bacteria", "eukaryota", "virus"))
WB.sub <- arrange(WB.sub, taxon_type)
WB.sub$taxon_name = factor(WB.sub$taxon_name, levels = unique(WB.sub$taxon_name))
unique(WB.sub$Primary.Diagnosis)
unique(WB.sub$Admission.Type)
WB.sub$Admission.Type[is.na(WB.sub$Admission.Type)] <- "not admitted"
WB.sub$Admission.Type = factor(WB.sub$Admission.Type, levels = c("not admitted", "fever"))
WB.sub <- arrange(WB.sub,  Admission.Type)
head(WB.sub)
WB.sub$sampleID <- NA
for (i in 1:length(unique(WB.sub$sample_name[WB.sub$Admission.Type=="not admitted" & WB.sub$Water.Control==" No"]))){
  WB.sub$sampleID[WB.sub$sample_name==unique(WB.sub$sample_name[WB.sub$Admission.Type=="not admitted" & WB.sub$Water.Control==" No"])[i]] <- paste0("control-healthy-", i)
}

for (i in 1:length(unique(WB.sub$sample_name[WB.sub$Admission.Type=="not admitted" & WB.sub$Water.Control=="Yes"]))){
  WB.sub$sampleID[WB.sub$sample_name==unique(WB.sub$sample_name[WB.sub$Admission.Type=="not admitted" & WB.sub$Water.Control=="Yes"])[i]] <- paste0("control-water-", i)
}


for (i in 1:length(unique(WB.sub$sample_name[WB.sub$Admission.Type=="fever"]))){
  WB.sub$sampleID[WB.sub$sample_name==unique(WB.sub$sample_name[WB.sub$Admission.Type=="fever"])[i]] <- paste0("fever-", i)
}

#WB.sub <- arrange(WB.sub, sampleID)
WB.sub$sampleID = factor(WB.sub$sampleID, levels=unique(WB.sub$sampleID))

p6 <- ggplot(data=WB.sub) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + theme_bw() + 
  theme(axis.title = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=90, size=4.5), panel.grid = element_blank(), 
        axis.text.y=element_text(size=3.3), strip.background = element_blank(), panel.spacing=unit(0, "lines"), legend.text = element_text(size=5)) +
  scale_fill_gradient(low ="yellow",  high = "red") + scale_x_discrete(position = "top") 
print(p6)

ggsave(file = "GCE_WB_all_pathogens.pdf",
       units="mm",  
       width=70, 
       height=150, 
       scale=3, 
       dpi=300)


SP.sub <- subset(heat.sub, Sample.Type=="serum" | Sample.Type=="plasma" | Sample.Type=="water")
SP.sub$taxon_type = factor(SP.sub$taxon_type, levels = c("bacteria", "eukaryota", "virus"))
SP.sub <- arrange(SP.sub, taxon_type)
SP.sub$taxon_name = factor(SP.sub$taxon_name, levels = unique(SP.sub$taxon_name))
unique(SP.sub$Primary.Diagnosis)
unique(SP.sub$Admission.Type)
unique(SP.sub$Study.Name)
SP.sub$Admission.Type[is.na(SP.sub$Admission.Type)] <- "not admitted"
SP.sub$Admission.Type = factor(SP.sub$Admission.Type, levels = c("not admitted", "fever"))
SP.sub <- arrange(SP.sub,  Admission.Type)
head(SP.sub)

SP.sub$Study.Name = factor(SP.sub$Study.Name, levels=c("water", "CZB", "ZORA", "Arbovirus"))
SP.sub <- arrange(SP.sub, Study.Name)
SP.sub$sampleID <- NA
for (i in 1:length(unique(SP.sub$sample_name[SP.sub$Study.Name=="ZORA"]))){
  SP.sub$sampleID[SP.sub$sample_name==unique(SP.sub$sample_name[SP.sub$Study.Name=="ZORA"])[i]] <- paste0("control-healthy-", i)
}

for (i in 1:length(unique(SP.sub$sample_name[SP.sub$Admission.Type=="not admitted" & SP.sub$Water.Control=="Yes"]))){
  SP.sub$sampleID[SP.sub$sample_name==unique(SP.sub$sample_name[SP.sub$Admission.Type=="not admitted" & SP.sub$Water.Control=="Yes"])[i]] <- paste0("control-water-", i)
}


for (i in 1:length(unique(SP.sub$sample_name[SP.sub$Admission.Type=="fever"]))){
  SP.sub$sampleID[SP.sub$sample_name==unique(SP.sub$sample_name[SP.sub$Admission.Type=="fever"])[i]] <- paste0("fever-", i)
}


SP.sub$sampleID = factor(SP.sub$sampleID, levels=unique(SP.sub$sampleID))

p7 <- ggplot(data=SP.sub) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + theme_bw() + 
  theme(axis.title = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=90, size=4.5), panel.grid = element_blank(), 
        axis.text.y=element_text(size=3.3), strip.background = element_blank(), panel.spacing=unit(0, "lines"), legend.text = element_text(size=5)) +
  scale_fill_gradient(low ="yellow",  high = "red") + scale_x_discrete(position = "top") 
print(p7)

ggsave(file = "GCE_SP_all_pathogens.pdf",
       units="mm",  
       width=70, 
       height=150, 
       scale=3, 
       dpi=300)



#now look at prevalence of key pathogens - only in the NP swabs
head(NP.sub)
NP.sum <- ddply(NP.sub, .(taxon_type, taxon_name), summarize, N_pos = length(sample_name))
NP.sum <- arrange(NP.sum, desc(N_pos))
NP.sum$N <- length(unique(NP.sub$sample_name))
NP.sum$prevalence <- NP.sum$N_pos/NP.sum$N
head(NP.sum)
NP.sum <- arrange(NP.sum, taxon_type, desc(prevalence))
NP.sum$taxon_name = factor(NP.sum$taxon_name, levels=unique(NP.sum$taxon_name))


p8 <- ggplot(data=NP.sum) + geom_bar(aes(x=taxon_name, y=prevalence, fill=taxon_type), stat = "identity") + theme_bw() + coord_flip() + #scale_y_discrete(position = "top") +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text( size=3), legend.position = c(.85,.85), legend.title = element_blank()) 
print(p8)

ggsave(file = "GCE_NP_prevalence.pdf",
       units="mm",  
       width=50, 
       height=150, 
       scale=3, 
       dpi=300)


#look at age-prevalence of viruses: Orthopneumovirus, Enterovirus, Cytomegalovirus
vir.sum = subset(NP.sum, taxon_type=="virus")
vir.sub = subset(NP.sub, taxon_type=="virus")
vir.sub = droplevels(vir.sub)
mean.vir.age <- ddply(vir.sub, .(taxon_name), summarize, mean.age= mean(Host.Age))
mean.vir.age <- arrange(mean.vir.age, desc(mean.age))

vir.sub$taxon_name=factor(vir.sub$taxon_name, levels= unique(mean.vir.age$taxon_name) )

#get age distribution of each infection
beeswarm(Host.Age~taxon_name, data=vir.sub)

head(vir.sub)
age.dist <- ddply(vir.sub, .(Host.Age), summarize, Ntot=length(sample_name))
vir.age <- ddply(vir.sub, .(Host.Age, taxon_name), summarize, Npos=length(sample_name))
vir.age$Ntot = NA
for (i in 1:length(age.dist$Host.Age)){
  vir.age$Ntot[vir.age$Host.Age==age.dist$Host.Age[i]]<- age.dist$Ntot[i]
}
vir.age$prevalence= vir.age$Npos/vir.age$Ntot
vir.age.plot = subset(vir.age, taxon_name=="Orthopneumovirus" | taxon_name== "Enterovirus" )#| taxon_name== "Cytomegalovirus")
p9 <- ggplot(data=vir.age.plot) + geom_point(aes(x=Host.Age, y=prevalence, size = Ntot, color=taxon_name), pch=1) + geom_line(aes(x=Host.Age, y=prevalence, color=taxon_name)) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank()) + xlab("host age") + facet_grid(taxon_name~.) + coord_cartesian(expand=F)
print(p9)

#why do enteroviruses show up in the respiratory tract?

#ggsave(file = "GCE_NP_age-prevalence.pdf",
#      units="mm",  
#     width=50, 
#    height=150, 
#   scale=3, 
#  dpi=300)

#now compare the positive controls
head(heat.sub)
unique(heat.sub$Primary.Diagnosis)
pos.control <- subset(heat.sub, !is.na(Primary.Diagnosis))
pos.control <- subset(pos.control, Primary.Diagnosis!="negative" & Primary.Diagnosis!="healthy")
unique(pos.control$sample_name)
unique(pos.control$Primary.Diagnosis)
unique(pos.control$Sample.Type)
pos.control$Sample.Type <- factor(pos.control$Sample.Type, levels=c("nasopharyngeal swab","whole blood", "serum" ))
pos.control <- arrange(pos.control, Sample.Type)
pos.control$Sample.Type <- as.character(pos.control$Sample.Type)
pos.control$sampleID <- NA
for (i in 1:length(pos.control$sampleID[pos.control$Sample.Type=="nasopharyngeal swab"])){
  pos.control$sampleID[pos.control$Sample.Type=="nasopharyngeal swab"][i] <- paste0("NP-swab-", seq(1, length(pos.control$sampleID[pos.control$Sample.Type=="nasopharyngeal swab"]), 1))[i]  
}
for (i in 1:length( pos.control$sampleID[pos.control$Sample.Type=="whole blood"])){
  pos.control$sampleID[pos.control$Sample.Type=="whole blood"][i] <- paste0("whole-blood-", seq(1, length( pos.control$sampleID[pos.control$Sample.Type=="whole blood"]), 1))[i]  
}
for (i in 1:length(pos.control$sampleID[pos.control$Sample.Type=="serum"])){
  pos.control$sampleID[pos.control$Sample.Type=="serum"][i] <- paste0("serum-", seq(1, length(pos.control$sampleID[pos.control$Sample.Type=="serum"]), 1))[i]  
}



p.pos.control <- ggplot(data=pos.control) + geom_tile(aes(x=sampleID, y= taxon_name, fill=log(NT_rpm))) + theme_bw() + 
  theme(axis.title = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=90, size=4.5), panel.grid = element_blank(), 
        axis.text=element_text(size=3.3), strip.background = element_blank(), panel.spacing=unit(0, "lines"), legend.text = element_text(size=5)) +
  scale_fill_gradient(low ="yellow",  high = "red") + scale_x_discrete(position = "top") 
print(p.pos.control)

ggsave(file = "GCE_pos.control.pdf",
       units="mm",  
       width=50, 
       height=150, 
       scale=3, 
       dpi=300)

#and load summary
pos.dat = read.csv(file = "pos.control.correct.csv", header=T, stringsAsFactors = F)
head(pos.dat)

#and plot
pos.sum <- ddply(pos.dat, .(sample.type, pathogen_type), summarise, Npos=sum(correct_diagnosis), Ntot=length(correct_diagnosis))
pos.sum$prop <- pos.sum$Npos/pos.sum$Ntot
pos.sum$pathogen_type[pos.sum$pathogen_type=="protozoan"] <- "eukaryota"
pos.sum$sample.type <- as.character(pos.sum$sample.type)
pos.sum$sample.type[pos.sum$sample.type=="nasopharyngeal swab"] <- "nasopharyngeal\nswab"
pos.sum$sample.type[pos.sum$sample.type=="whole blood" & pos.sum$pathogen_type=="bacteria"] <- "whole blood\n(bacteria)"
pos.sum$sample.type[pos.sum$sample.type=="whole blood" & pos.sum$pathogen_type=="eukaryota"] <- "whole blood\n(eukaryota)"
pos.sum$sample.type[pos.sum$sample.type=="whole blood" & pos.sum$pathogen_type=="virus"] <- "whole blood\n(virus)"
pos.sum$sample.type = factor(pos.sum$sample.type, levels = c("nasopharyngeal\nswab", "whole blood\n(bacteria)", "whole blood\n(eukaryota)", "whole blood\n(virus)", "serum"))
pos.sum$N_label = paste0("N=", pos.sum$Ntot)

pos.correct <- ggplot(data=pos.sum) + geom_bar(aes(x=sample.type, y=prop, fill=pathogen_type), stat="identity", position="dodge") +ylab("proportion diagnosed correctly") + coord_cartesian(ylim=c(0,1))+
  theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(),panel.grid = element_blank(), legend.position = c(.75,.95), legend.text = element_text(size=12), axis.title = element_text(size=16), axis.text = element_text(size = 12), legend.direction = "horizontal") +
  geom_label(aes(x=sample.type, y=.92, label=N_label), size=3)
print(pos.correct)
ggsave(file = "GCE_pos.control_comp.pdf",
       units="mm",  
       width=90, 
       height=50, 
       scale=3, 
       dpi=300)

#for poster
pos.correct.poster <- ggplot(data=pos.sum) + geom_bar(aes(x=sample.type, y=prop, fill=pathogen_type), stat="identity", position="dodge") +ylab("proportion diagnosed correctly") + coord_cartesian(ylim=c(0,1))+
  theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(),panel.grid = element_blank(), legend.position ="top", legend.text = element_text(size=16), axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 16), legend.direction = "horizontal") +
  geom_label(aes(x=sample.type, y=.92, label=N_label), size=4)
print(pos.correct.poster)
ggsave(file = "GCE_pos.control_poster.pdf",
       units="mm",  
       width=80, 
       height=50, 
       scale=3, 
       dpi=300)


#and prop non-host reads per microbe type
head(heat.sub)
read.sub <- ddply(heat.sub, .(Sample.Type, sample_name, taxon_type), summarize, total_reads = sum(NT_r))

#and total reads per sample
read.tot <- ddply(read.sub, .(Sample.Type, sample_name), summarize, N_total_reads = sum(total_reads))
head(read.tot)

read.sub$N_total_reads <-NA
for(i in 1:length(read.tot$Sample.Type)){
  read.sub$N_total_reads[read.sub$sample_name==read.tot$sample_name[i]] <- read.tot$N_total_reads[i]
}
read.sub$prop <- read.sub$total_reads/read.sub$N_total_reads
head(read.sub)

#sort
read.sub$Sample.Type = factor(read.sub$Sample.Type, levels=c("water", "whole blood", "plasma", "serum", "nasopharyngeal swab"))
read.sub <- arrange(read.sub, Sample.Type)
read.sub$sample_name = factor(read.sub$sample_name, levels = unique(read.sub$sample_name))


p.out <- ggplot(data=read.sub) + geom_bar(aes(x=sample_name, y= prop, fill=taxon_type, color=Sample.Type), stat = "identity", position = "stack")+theme_bw()+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=18), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        legend.position = "top", legend.title=element_blank(), legend.direction = "horizontal") + guides(fill = guide_legend(nrow = 1))
print(p.out)

ggsave(file = "GCE_prop_reads_taxon_type_spacer.pdf",
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)


p.out1 <- ggplot(data=read.sub) + geom_bar(aes(x=sample_name, y= prop, fill=taxon_type), stat = "identity", position = "stack")+theme_bw()+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=18), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        legend.position = "top", legend.title=element_blank(), legend.direction = "horizontal") + guides(fill = guide_legend(nrow = 1))
print(p.out1)

ggsave(file = "GCE_prop_reads_taxon_type.pdf",
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)



p.test <- ggplot(data=read.sub) + geom_bar(aes(x=sample_name, y= prop, fill=taxon_type), stat = "identity", position = "stack")+theme_bw()+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=18), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x = element_text(size=2, angle=90),
        legend.position = "top", legend.title=element_blank(), legend.direction = "horizontal") + guides(fill = guide_legend(nrow = 1))
print(p.test)
ggsave(file = "GCE_prop_reads_with_name.pdf",
       units="mm",  
       width=140, 
       height=50, 
       scale=3, 
       dpi=300)

