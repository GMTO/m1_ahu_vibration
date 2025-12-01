# Repository with scripts to calculate the M1 AHU transfer functions

The script `calc_m1ahu_tfs.mlx` combines a state-space model representing the GMT's structural dynamics, the mount axes feedback controller transfer functions, optical sensitivity matrices, and rejection transfer functions to obtain the frequency responses between AHU force (and torque) inputs and segment tip-tilt and piston, and corresponding wavefront error at the exit pupil. The goal is to provide a transfer function matrix representing the telescope's sensitivity to forces applied to nodes at the M1 AHU (air handling units) locations.
The script to calculate the M1 AHU frequency responses relies on the following data files:
- `odc_fdr2023_z30_o`: It provides a struct object "o" with mount subsystem models, including the mount controller transfer functions.
- `lom_tt_dt` and `D_seg_piston_dt`: Those files provide linear approximations of the transformations from the M1 and M2 rigid-body motions to segment tip-tilt and piston at the telescope exit pupil.
- `modal_state_space_model_2ndOrder`: The file contains the telescope structural dynamics model represented in the second-order modal form. The variable `modelDescription` provides a detailed description of the model.

After running the script, the state-space model matrices are saved as `G_20251110_1617z30zeta<.>pct_ahu_ss.mat`, and the output data `G_20251110_1617z30zeta<.>pct_ahu_fr.mat` file provides the frequency response model for each observing mode: `H_NS_ttp`, `H_NGAO_ttp`, and `H_LTAO_ttp`.
