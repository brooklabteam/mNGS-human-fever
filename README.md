# mNGS-human-Fever

This repo walks through some simple analyses of the CZID sequencing data from the undiagnosed human fever project uploaded [here](https://czid.org/my_data?currentDisplay=table&currentTab=samples&mapSidebarTab=summary&projectId=516&showFilters=true&showStats=true&updatedAt=2022-01-06T04%3A24%3A21.974Z&workflow=short-read-mngs). You can start by playing around on CZID to get a feel for the data that are available. Then, you can move analyses over here into R.

First, on CZID, try highlighting all the samples, and in the upper righthand corner of the sample list, click 'Download' - you will see several options for types of files to download. First, try clicking the 'Samples Overview' csv file for download -- we have stored that here in the 'data' subfolder under the name "gce_sample_summary.csv". In the "R-code" subfolder, you will find a script "process_plot_summary.R" that walks through how to visualize this output.

There is also the option to download "Sample Metadata", which is user-uploaded metadata that has been uploaded with the sequences. Depending on who did the upload, however, these data are often incomplete. Here, we will supply our own metadata file instead. 

Next, if you want to download the pathogen hits associated with each sample, try selecting all the samples, then building a 'heatmap', then download the corresponding csv file from the heatmap.
