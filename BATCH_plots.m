function BATCH_plots(ModelName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These plots can only be generated after the forward modeling has been performed and the run subdirectories (e.g., 
% "Archive-Traplike_1.00Ref" populated with run output
% This is most easily accomplished using the included "BATCH_run" script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model names are character strings derived from run data files, for example:

% File name                   Model name

% Trapridge_1.00Ref.dat       Trapridge_1.00Ref
% Traplike_1.00UR_ML_LR.dat   Traplike_1.00UR_ML_LR

% Plot results are stored subdirectory "PLOTS" in run archive subdirectories such as "Archive-Traplike_1.00Ref"

% Plot 1 : S_0 surface stratification simulation results

fprintf(1,'\n=>> main_plot_surface_stratification_S0(%s)\n', ModelName);

main_plot_surface_stratification_S0(ModelName)

% Plot 2 : S_1 strain foliation simulation results

fprintf(1,'\n=>> main_calculate_and_stereoplot_strain_foliation_S1(%s)\n', ModelName);

main_calculate_and_stereoplot_strain_foliation_S1(ModelName)

fprintf(1,'\n=>> main_calculate_and_plot_strain_foliation_S1_map(%s)\n', ModelName);

main_calculate_and_plot_strain_foliation_S1_map(ModelName)

% Plot 3: Folding parameter profile 

fprintf(1,'\n=>> main_plot_folding_parameter_profile(%s)\n', ModelName);

main_plot_folding_parameter_profile(ModelName)
