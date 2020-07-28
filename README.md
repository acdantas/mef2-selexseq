"# mef2-selexseq" 
Scripts for the analyses of MEF2 data.

1) For analysis of raw data, selex package is used as described in the manuscript

2) Once kmer table of affinities is obtained, scripts were ran as follow to analyze and plot data:

set-parameters.R

get-kmer-analysis-1.R

get-kmer-analysis-2.R

3) Filtered data generated (./data/mlrseq) were used to run MLR according to an approach previously published and described on the manuscript. 

4) After obtaining MLR results (sumamry.txt), the results were plotted with:

get-MLR-analysis.R

