Included in this repository are raw data files and code related to the paper "The Density of Braun's Lipoprotein Determines Vesicle Production in E. coli"

**Lpp_Vesicles_data_fig3**.xlsx: Includes measurements of vesicles/cell used to determine y coordinates in panel b of figure 3.

**Lpp_XY_figS6_S7.mat**: Lpp density and vesicles/cell values used for data fitting in figs S6 and S7

**Lpp_antiCluster_figS2.py**: Code used to generate crosslink locations and circle sizes. 
This implements a minimum separation distance for crosslinks, shown in fig S2.

**Lpp_cluster_fig2_figS2.py**: Code used to generate crosslink locations and circle sizes. 
This implements a clustering scheme for crosslinks, shown in fig S2.
Setting cluster size to 0 prevents clusters from forming, and produced results shown in Fig 2.

**Lpp_scaling_fit_fig_S6_S7.m**: Takes data from Lpp_XY... and runs RC_KP_Lpp... and fits for the parameters B and D,
which are shown in fig S6, fig S7, and text S2.

**RC_KP_Lpp_fig_S6_S7.m**: Generates Rc values for a given choice of parameters B or D. Used for data fitting.

**Raw_qPCR.xlsx**: Raw qPCR data

**pTetGFP_Stationary2_figS5.xlsx**: Fluorescence data from pTetGFP strain taken in stationary phase.

**pTetGFP_Stationary_figS5.xlsx**: Fluorescence data from pTetGFP strain taken in stationary phase.

**pTetGFP_Exponential_Fig3.xlsx**: Fluorescence data from pTetGFP strain taken in exponential phase.

**qPCR_Process_Lpp**: Reads data from Raw_qPCR.xlsx and calculates relative density of Lpp at each induction level.
