# GlacierStructureModeling
Simulation of meso-scale structures in surge-type glaciers

This repository contains MATLAB scripts and data files used in generating figures and modeling results in two papers that we hope to publish in the Journal of Geophysical Research - Earth Surface

Hambrey, M.J., & Clarke, G.K.C. 2019. Structural evolution during cyclic glacier surges: 1. Structural glaciology of Trapridge
  Glacier, Yukon, Canada. Journal of Geophysical Research--Earth Surface.
  
Clarke, G.K.C., & Hambrey, M.J. 2019. Structural evolution during cyclic glacier surges: 2. Numerical modeling. Journal of
  Geophysical Research--Earth Surface.
  
No plots or run output are contained in this repository so it is necessary to run the various programs to generate output and plot files for each model considered. Depending on the model selected, the archived output for a single model run is roughly 250 GB and the run times will exceed 24 hours on a conventional workstation.

The subdirectories  ./Archive-Traplike_1.00Ref,  ./Archive-Trapridge_1.00Ref, etc
refer to models 'Traplike_1.00Ref', 'Trapridge_1.00Ref', etc. and contain data files relevant to the particular model. For example, the mass balance forcing (e.g., b_dot_v2016-09-10.mat), sliding mask (e.g., sliding_mask_and_surge_axis_v2016-09-10.mat), bed surface digital elevation model (e.g.,Traplike_BS_dem_v2016-09-09.mat). The file SiteList.dat provides information about measurement sites (actual or imagined, including sites within the body of the glacier). The file SourceList.dat gives information about the location of surface inputs for rock debris (this is used in the medial moraine model). The file PlotList.dat contains information used by the plotting scripts (e.g., which sites to plot and whether hardcopy PDF output is desired).

Within the top directory of the repository (i.e., the same directory as for this README.md file) are data files for each of the models, for example, Traplike_1.00Ref.dat, Traplike_1.00NoSurge.dat, Trapridge_1.00Ref, etc.

Getting started:

You can run any given model using the BATCH_run.dat file which calls a series of "main" programs starting with "main" which runs the thermomechanical ice dynamics file and stores run snapshots at 0.01 year intervals for a single surge cycle. These snapshots are stored in the subdirectories ./Snapshots with the ./Archive-xxx   subdirectories. Next bottom-up and top-down particle tracking are performed, etc.

To launch a batch run select the desired model, e.g.,'Traplike_1.00Ref'  and enter BATCH_run('Traplike_1.00Ref') at the MATLAB prompt. I usually choose to do these runs as background jobs so in Linux here is how one might proceed:

Inside a terminal:        

nice -n 20 matlab -nodesktop
BATCH_run('Traplike_1.00Ref')
                         
Note the BATCH_run will take several days to complete. First the forward model is run. Following completion of the forward modeling run: particle tracking and medial moraine generation are computationally intensive 
                        
Results of bottom-up particle tracking are stored in the file UP_TRACKS.mat within the Archive-xxx subdirectory
Results of top-down particle tracking are stored in the file DOWN_TRACKS.mat in the same Archive-xxx subdirectory
                         
Once the BATCH_run is completed you can BATCH process some plots using

BATCH_plot('Traplike_1.00Ref')
                         
