# M1 AHU Vibration Analysis Scripts

This repository provides scripts developed to assess M1 AHU (air handling unit) vibration data and its effect on the telescope performance.

The main scripts are:
- `plot_lc_locations.mlx`: plots the load cell locations used to measure the AHU vibration forces.
- `proc_ahu_lc_dt.mlx`: plots the AHU force data and calculate the vibration force power spectrum density (PSD). Besides plotting and reporting relevant information regarding the vibration data, the last cell of the script enables the user to save the PSD data so that it can be used in further analysis (see `m1_ahu_vib_analysis.mlx`)
- `compute_m1fan_rtf.mlx`: computes the frequency response of the telescope transfer function from the forces applied at the M1 fan locations to the segment tipt-tilt and piston at the exit pupil. Those data are crucial to estimate the impact of the vibration forces in the telescope performance.
- `m1_ahu_vib_analysis.mlx`: provides an estimate of the wavefront error induced by the M1 AHU vibration forces. The script employs the frequency response of the telescope transfer function matrix (`compute_m1fan_rtf.mlx`) and the PSD of AHU vibration forces. Pre-computed frequency response files and AHU vibration PSD data are available on DROBO (`smb://Drobo-IM.gmto.org/im/Integrated Modeling/IM-studies/m1_ahu_vib/`).
- `m1_ata_vib_analysis.mlx`: provides tip-tilt and piston responses to white-noise inputs applied to the M1 AHU nodes. The script employs pre-computed frequency response of the telescope transfer function matrix (`m1sXfanTF_4ataPSD.mat`). Those pre-computed frequency response files are also available on DROBO (`smb://Drobo-IM.gmto.org/im/Integrated Modeling/IM-studies/m1_ahu_vib/`). At last, the script computes the wavefront error impact due to the forces computed by ATA in 2018. The PSD of the AHU vibration data computed by ATA is available from file `Papst1_10.csv` saved at the same folder.
