function BATCH_run(ModelName)

% Model names are character strings derived from run data files such as

% File name                   Model name

% Trapridge_1.00Ref.dat       Trapridge_1.00Ref
% Traplike_1.00UR_ML_LR.dat   Traplike_1.00UR_ML_LR

% Results are stored in run archive subdirectories such as "Archive-Traplike_1.00Ref"

% Step 1: Run the thermodynamic flow model

main(ModelName);

% Step 2: Particle tracking from downstream site to upstream source and back

main_bottom_up_particle_tracking(ModelName);
main_top_down_particle_tracking_etc(ModelName);
main_track_closure_check(ModelName);

% Step 3: Postprocessing of down-tracking results

% Gradient of velocity gradient, gradient of displacement gradient and calculation of rank 3 invariants
% that characterize folding

main_calculate_and_save_grad_L(ModelName);
main_calculate_and_save_F_par_structure(ModelName);
main_calculate_folding_parameters(ModelName);

% Step 4: Crack analysis 

main_calculate_and_stereoplot_crack_data(ModelName);
main_calculate_tabulate_and_plot_crack_density(ModelName);

% Step 5: Moraine track calculations

main_construct_moraine_debris_emergence_archive(ModelName);
main_calculate_moraine_tracks(ModelName);


