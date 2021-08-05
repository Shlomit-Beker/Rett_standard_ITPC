# Rett_standard_ITPC
codes for analyzing ITPC in Tufi et al Rett paper
The codes here are for the analysis of ITPC data in two different ways: 1. for each frequency in the range, collapsed across all time points 2. for each frequency and time points. 
the coes in the order of execution are: 
1. Main_Rett: creates the data structure that is used for the analysis, from the individual subject epoched structures.
2. PL_rett: either run after Main_Rett or upload AllData. Calculates ITPC across time points. 
3. timeFreq_Rett: either run after Main_Rett or upload AllData. Calculates ITPC frequency-time ITPC maps
4. it_plotSmoothTFR_rett: plot any ITPC data on 3d image. 

Other functions: 
1. All_Data_ds: downsampling the AllData structure, for the time-frequency ITPC.
2. FDR_clusBased_permuTest: runs cluster-based permutation test and plot in transparent/bold colors for pre-post correction, respectively. 
3. basewave4: compute the time-frequency representation with Morlet wavelet (written/editted by Peter Lakatos). 
