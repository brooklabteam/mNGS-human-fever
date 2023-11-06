# mNGS-human-Fever

This repo walks through some simple analyses of the CZID sequencing data from the undiagnosed human fever project uploaded [here](https://czid.org/my_data?currentDisplay=table&currentTab=samples&mapSidebarTab=summary&projectId=516&showFilters=true&showStats=true&updatedAt=2022-01-06T04%3A24%3A21.974Z&workflow=short-read-mngs). You can start by playing around on CZID to get a feel for the data that are available. Then, you can move analyses over here into R.

First, on CZID, try highlighting all the samples, and in the upper righthand corner of the sample list, click 'Download' - you will see several options for types of files to download. First, try clicking the 'Samples Overview' csv file for download -- we have stored that here in the 'data' subfolder under the name "gce_sample_summary.csv". 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/download-types.png  width="500" height="500">


In the "R-code" subfolder, you will find a script "process_plot_summary.R" that walks through how to visualize this output. The resulting plots start with "QC_" and can be found in the "figures" folder.

On CZID, there is also the option to download "Sample Metadata", which is user-uploaded metadata that has been uploaded with the sequences. Depending on who did the upload, however, these data are often incomplete. Here, we will supply our own metadata file instead - you can find this one in the "data" subfolder under name "GCE-human-metadat.csv". 

Also, on CZID, if you want to download the pathogen hits associated with each sample, try selecting all the samples, then building a 'heatmap', then download the corresponding csv file from the heatmap. You have the option to download the heatmap after filtering or before. There are a few subsets of these heatmaps that are referenced in the 'GCE_human_analyses.R' script. I downloaded them as separate subsets to help with the addition of higher level taxonomic groupings (CZID only downloads pathogen name for example), as well as sample type.

Finally, it is possible to download genomic data for a single sample (or several samples simultaneously), as well as the heatmap metrics. Take a look at [this sample](https://czid.org/samples/23806) as an example. If you click on "Metapneumovirus", you will see the converage visualization and the two contigs that were constructed. You can download those directly as .fasta files by clicking on the cloud icon to the right. 

<img src=https://github.com/brooklabteam/mNGS-human-fever/blob/main/guide-pics/coverage-visualization.png  width="700" height="300">



