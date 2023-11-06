# mNGS-human-Fever

This repo walks through some simple analyses of the CZID sequencing data from the undiagnosed human fever project uploaded [here](https://czid.org/my_data?currentDisplay=table&currentTab=samples&mapSidebarTab=summary&projectId=516&showFilters=true&showStats=true&updatedAt=2022-01-06T04%3A24%3A21.974Z&workflow=short-read-mngs). You can start by playing around on CZID to get a feel for the data that are available. Then, you can move analyses over here into R.

First, on CZID, try highlighting all the samples, and in the upper righthand corner of the sample list, click 'Download' - you will see several options for types of files to download. First, try clicking the 'Samples Overview' csv file for download -- we have stored that here in the 'data' subfolder under the name "gce_sample_summary.csv". 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/download-types.png  height="300">


In the "R-code" subfolder, you will find a script "[process_plot_summary.R](https://github.com/brooklabteam/mNGS-human-fever/blob/main/R-code/process_plot_summary.R)" that walks through how to visualize this output. The resulting plots start with "QC_" and can be found in the "figures" folder.

Here are a few examples of what is produced - proportion of reads cleared by each filtration step in CZID:
<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/figures/QC_GCE_prop_reads_spacer.png  width="500">

QC metrics by sample type:
<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/figures/QC_GCE_by_sample_type.png  width="300">


On CZID, there is also the option to download "Sample Metadata", which is user-uploaded metadata that has been uploaded with the sequences. Depending on who did the upload, however, these data are often incomplete. Here, we will supply our own metadata file instead - you can find this one in the "data" subfolder under name "GCE-human-metadat.csv". 

Also, on CZID, if you want to download the pathogen hits associated with each sample, try selecting all the samples, then building a 'heatmap', then download the corresponding csv file from the heatmap. You have the option to download the heatmap after filtering or before. There are a few subsets of these heatmaps that are referenced in the '[GCE_human_analyses.R](https://github.com/brooklabteam/mNGS-human-fever/blob/main/R-code/GCE_human_analyses.R)' script. I downloaded them as separate subsets to help with the addition of higher level taxonomic groupings (CZID only downloads pathogen name for example), as well as sample type.

Among other things, these show you how to make heatmaps like the following:
<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/figures/GCE_NP_swab_all_pathogens_heatmap.png  width="300">

And how to summarise the data like so:

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/figures/GCE_prop_reads_taxon_type.png  width="500">


Finally, it is possible to download genomic data for a single sample (or several samples simultaneously), as well as the heatmap metrics. Take a look at [this sample](https://czid.org/samples/23806) as an example. If you click on "Metapneumovirus", you will see the coverage visualization plot below and the two contigs that were constructed. You can download those directly as .fasta files by clicking on the cloud icon to the right. 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/coverage-visualization.png  width="600">

You can also download all of the non-host contigs from this sample, (meaning all of the contigs that were assembled by de novo assembly in CZID after all of the filtration steps) by clicking on the "Download" button with the cloud in the upper right of the screen. 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/all-contig-download.png  width="550">

These contigs can then be BLASTed (either manually or via a command line script) to NCBI. Occassionally, CZID assembles a contig correctly but gets the BLAST link wrong, so a manual BLAST to a curated reference database can be helpful. Gwen used this project to identify contigs that were hits to bat CoVs her [paper](https://www.frontiersin.org/articles/10.3389/fpubh.2022.786060/full), following this pipeline [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/contig-blast-directions.md).

Additionally, you can download the unmapped non-host reads (meaning all reads that passed filter but not yet assembled into contigs) by clicking on the "View Results Folder" line from the "Download" button and scrolling to the very bottom. This two nonhost fastq files are Illumina paired end reads. 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/download-raw-reads.png  width="400">


You can use these to do your own de novo assembly or to map reads to a contig to see the coverage depth (e.g. support) for a particular genome (essentially recreating the coverage plot above). I walk through a couple examples of these paired covereage plots using Christian's script, comparing coverage from [MSSPE](https://www.nature.com/articles/s41564-019-0637-9) and mNGS in the script '' located in the R-scripts folder. A few plots found in the "figures" sub-folder are then produced. This script compares the HCoV-HKU1 coverage from mNGS (see [here](https://czid.org/samples/23804) and click on the 'Consensus Genome' tab) vs. that with MSSPE (see [here](https://czid.org/samples/358165)).

This highlights one other feature of CZID--you can build a 'consensus genome' where it maps all raw reads to the closet hit in GenBank and tries to build a consensus genome. You will see this as an option when you upload new samples. You can download the data needed to produce above by clicking the 'Download All' button (top right with cloud) on each sample page in CZID. 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/download-all-vis.png  width="400">


Clicking on the 'Download All' option downloads a zipped file that actually has a read depth plot as well as the non-host reads (fq.gz files) and several output files. The 'samtools_depth.txt' file that is produced gives you the read depth across the entire genome. I included these two samtools files and teh two .tsv report files for mngs and msspe in the 'data' subfolder -- see '' script to process these into a coverage plot like this:


<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/figures/mNGS_MSSPE_comparison_read_depth.png  width="500">

The interesting thing about above is that, by just looking at the raw read depth plot, it looks like MSSPE is doing great. But when you actually plot reads per million, we see that that mNGS actuallyhas higher reads per million in many places in the genome, but that there were just A LOT of reads produced in the MSSPE run (this makes sense as it was run on a NovaSeq while the mNGS was run on a NextSeq). However, we can say that the coverage has greater breadth in the case of MSSPE, meaning no dropouts with no resolution across the genome. This is a value that us reported in the .tsv report files as "Coverage >= 1x (%)" -- you'll see it is 99% for the MSSPE (meaning every basepair has support at 1X or higher), while for mNGS it is only 1.71%.


Finally, once you have genomes in hand, the Brook lab has lots of great resources for how to build Bayesian timetrees (e.g. [here](https://github.com/brooklabteam/cambodia-dengue-national)) or maximum likelihood phylogenies (e.g. [here](https://github.com/brooklabteam/cambodia-dengue-national/blob/main/figure-development/FigS18/Prep-ML-Tree.md)). There are also good how-to scripts for these in the [Mada-Bat-CoV repo](https://github.com/brooklabteam/Mada-Bat-CoV).

It's also worth taking a glance through the [SARS-CoV-2 genome curation repo](https://github.com/brooklabteam/SC2-genome-curation) to get an idea of other secrets hidden in CZID! There are some scripts in here that automate the download of the samtools files to produce read coverage plots across many sample types.






