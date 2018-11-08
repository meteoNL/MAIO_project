The results of combining the dropbox files that we got for GPS-data at each site, after cropping them (removing off-site measurements) and resolving issues with e.g. separators and decimal delineators, is found in [sitename]".txt".

The data processing files for GPS-data are "read_data_2.py" and "read_data_3.py". The second version of this file is the most convenient one, as it just does the data processing and shows the most relevant graphs. The third version additionally produces giant graphs with awefully big titles which were used in the presentation to show the intermediate steps.
These files use the ".txt" files with location data. The processing steps are described in detail in the report. 
Using the "read_data_2.py", the [sitename]".csv" files (separated by ; and as decimal delineator a ,) were generated. These contain 96 hour average velocities. For S9, there was a special case, as Roderik sent processed velocities to us due to some missing data in the dropbox files. The values processed by us are featured in "S9 (old).csv" and the combination of two processed velocity records is found in "S9.csv". For different results, different choices out of the two were made, therefore this was convenient. 

The csv-files are then read by the figure generating files (finding relations of time trends and correlations). In addition to this, the IWS results are read by the same files (except seasonal). These were separately processed and the results are found in numpy array files, "iws"[season][sitename]".npy". All three files below use the csv-files to read the GPS-velocity data from. 
The files "cumbal"[sitename]".txt" are the cumulative balances of the sites. These are read in the "balancevelcorr+bal_analyis.py"-file. This file is used to combine annual mass balance observations (corresponding time axis in "melttimes.txt") with the seasonal and annual velocity records. Here, trends in mass balance and correlations with velocity are calculated for the time lags. Additionally, the significance of trends is analyzed in this file. 
A separate file is "seasonal.py", which basically throws all velocity records of a site in discretized time buckets of about 8 days and then averages these to produce a graph. 
Then we have "vel_time_trends.py" for analyzing the time trends in velocity. 

The "alllocations.xlsx"-file contains annual velocities per site as calculated in the balancevelcorr and mass balance, to compute weighted areal mean values of these variables. We couldn't finish this analysis, but yellow marked cells indicate values which deviate from earlier analyses.  


Last but not least, there are a couple of folders: 
The figures in the report are found in the folder figures, with subdirectories for older figures. 
The "supplement_VDWal2015"-folder contains the supplementary material of the paper  R.  S.  W.  van  de  Wal,  C.  J.  P.  P.  Smeets,  W.  Boot,  M.  Stoffelen,  R.  van  Kampen,  S.  H.  Doyle, F. Wilhelms, M. R. van den Broeke, C. H. Reijmer, J. Oerlemans, and A. Hubbard.  Self-regulation of ice flow varies across the ablation area in south-west greenland. The Cryosphere, 9(2):603–611, apr 2015.  doi:10.5194/tc-9-603-2015
The folder "presentatie" contains the powerpoint presentation given on November 8th in class. 
The other folders contain some intermediate result steps of the IWS data, most notably used for figure 3 of our report. These files are maybe not so well organized. 