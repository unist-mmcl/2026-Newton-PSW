# 2026-Newton-psw
This repository contains the code used to generate the simulation results presented in the figures of Kang et al., “Manipulating molecules using partial optical standing waves,” Newton (2026). 


## How to use

To run a custom trajectory simulation of CS₂ molecules, open the project's root folder in MATLAB.

- To run a simulation with a single parameter set, enter `run_all` in MATLAB. This loads the default parameters from `STN_constants.mat` and accepts an optional struct to override individual settings. [run_all.m](run_all.m) · [STN_constants.mat](STN_constants.mat)
- To run simulations with multiple parameter sets, edit `param_sets` and enter `run_batch` in MATLAB. Each parameter set overrides the default settings for a single call to `run_all`. [run_batch.m](run_batch.m)

To reproduce the simulations for specific figures, set MATLAB's current folder to the corresponding folder listed below and run the indicated command. Each figure uses the constants file in its own folder and saves the results in that folder's 'result' directory.

| Figure | MATLAB folder | Command |
| --- | --- | --- |
| Fig. 4 | `fig4/` | [`run_all_fig4`](fig4/run_all_fig4.m) |
| Fig. 5 | `fig5/` | [`run_batch_Fig5`](fig5/run_batch_Fig5.m) |
| Fig. 6 | `fig6/` | [`run_all_fig6`](fig6/run_all_fig6.m) |
| Figs. S1 and S2 | `figS1S2/` | [`run_all_figS1S2`](figS1S2/run_all_figS1S2.m) |
| Fig. S3 | `figS3/` | [`run_batch_figS3`](figS3/run_batch_figS3.m) |
| Fig. S4 | `figS4/` | [`run_batch_figS4`](figS4/run_batch_figS4.m) |

See [Figure-specific simulations](#figure-specific-simulations) for details on the quantum-state batch simulations used for Figs. 5, S3, and S4.

## Simulation workflow

[run_all.m](run_all.m) runs the CS₂ trajectory simulations for the intensity pairs specified in the current parameter set. It adds the project root directory and the 'sim' folder to the MATLAB search path.

| Step | What runs |
| --- | --- |
| 0. Configure the parallel pool | Initializes or reuses a parallel pool with eight workers, limiting each worker to a maximum of three computational threads. |
| 1. Load constants | Loads `STN_constants.mat`, sets the default particle count `MAX` to 10 million per intensity pair, then applies fields from the optional `override` struct. The `label` field sets the result filename suffix. |
| 2. Set image processing parameters | Sets the electron recoil parameters, velocity-to-pixel conversion, and x/y scaling factors for the Voigt point-spread function. |
| 3. Locate the quantum-state data | Resolves `jm_data` to a text file in `data/`. Its state probabilities are used to sample each particle's rotational state. |
| 4. Load the alignment table | Loads `linearly.mat` and flattens its `(J, abs(M))` dimensions into rows so workers can look up intensity-dependent alignment values. |
| 5. Compute intensities and divide the work | `I1_list` and `I2_list` specify the energy per pulse in mJ. The code converts these energies to joules and computes the peak intensities in W/m² using a transverse Gaussian beam profile with waists `w0y` and `w0z` and a temporal Gaussian pulse envelope with an intensity FWHM of `tau`. The factor `atten1` accounts for transmission losses through the experimental window, and `atten2` corrects the peak intensity for deviations of the measured transverse beam profile from an ideal Gaussian. Set either factor to `1` to omit its correction. `max_I` gives the peak intensity in W/m² for a 1 mJ input pulse after both attenuation corrections. The code multiplies each input pulse energy, expressed in mJ, by `max_I` to obtain its peak intensity. For each intensity pair, the code schedules sub-batches of 5,000 particles. |
| 6. Run the sub-batches in parallel | Each `parfor` task calls [run_parfor_task.m](sim/run_parfor_task.m). It first calls [init_particles.m](sim/init_particles.m) to sample states, initial positions, velocities, and laser jitter. It then selects the alignment row for each particle and calls [time_integrate.m](sim/time_integrate.m). The integrator calls [axay_opt_partial.m](sim/axay_opt_partial.m) at each time step to compute the optical dipole acceleration and update the particle motion. |
| 7. Reassemble velocities | Concatenates the sub-batch results for each intensity pair into one matrix with columns `vx`, `vy`, and `vz` in m/s. The matrices are stored in `Cell_velocity` in the order of the intensity lists. |
| 8. Process images (optional) | Calls [image_processing.m](sim/image_processing.m) when enabled. It adds electron recoil to the simulated velocities, bins them on a pixel grid, convolves each histogram with a 2-D Voigt point-spread function, and computes x/y projections. |
| 9. Save results | Writes `Cell_velocity`, `I1_list` and `I2_list` (pulse energies in mJ), and `constants` to `result/result_<label>.mat`. If no label is supplied, the filename uses the `jm_data` basename. |

### Image processing during or after a run

Step 8 is commented out by default. Uncomment the call to `image_processing` to compute `img_matrix`, `sumx`, and `sumy` during the simulation. To include those variables in the result MAT file, add them to the `save` call in step 9.

If you skip step 8, you can process the saved velocities later without rerunning the simulation. From the project root, load a CS₂ result and run [Image_code_voigt.m](Image_code_voigt.m). Replace the example path with your result file.

```matlab
load(fullfile('result', 'result_example.mat'), 'Cell_velocity', 'I2_list', 'constants');
MAX = size(Cell_velocity{1}, 1);
rep_std = constants.rep_std;
rep_vel = constants.rep_vel;
vel_pix = constants.vel_pix;
atten_factorx = constants.atten_factorx;
atten_factory = constants.atten_factory;
Image_code_voigt
```

The script produces `img_matrix`, `sumx`, and `sumy` in the MATLAB workspace.

### Batch runs

[run_batch.m](run_batch.m) builds a cell array of parameter structs named `param_sets`. For each entry it calls `run_all(param_sets{k})`. Fields such as `jm_data`, `I1_list`, `I2_list`, `tau`, and `MAX` replace the corresponding settings for that run. Each call starts from the defaults in step 1, then applies that entry's overrides. Set a distinct `label` for each entry when several runs use the same `jm_data`.

### Figure-specific simulations

Each figure folder contains its own constants file and MATLAB entry point. Set the MATLAB current folder to the corresponding directory and run the command shown below. The figure runners save their results under that folder's `result/` directory.

| Figure | Constants file | MATLAB command |
| --- | --- | --- |
| Fig. 4 | [fig4_constants.mat](fig4/fig4_constants.mat) | `run_all_fig4` |
| Fig. 5 | [Fig5_constants.mat](fig5/Fig5_constants.mat) | `run_batch_Fig5` |
| Fig. 6 | [fig6_constants.mat](fig6/fig6_constants.mat) | `run_all_fig6` |
| Figs. S1 and S2 | [figS1S2_constants.mat](figS1S2/figS1S2_constants.mat) | `run_all_figS1S2` |
| Fig. S3 | [figS3_constants.mat](figS3/figS3_constants.mat) | `run_batch_figS3` |
| Fig. S4 | [figS4_constants.mat](figS4/figS4_constants.mat) | `run_batch_figS4` |

The Fig. S3 and Fig. S4 batch scripts override `jm_data` with the matching `alpha_JM_*.txt` file for each of the nine states (`00`, `20`, `21`, `22`, `40`, `41`, `42`, `43`, and `44`). Their individual `run_all` functions use the default distribution `alpha_T1_J=2n,j,m.txt` stored in the corresponding constants file.

Fig. 5 simulates nitric oxide (NO) and uses a separate alignment table, [linearly_NO.mat](fig5/linearly_NO.mat), built from NO alignment data. [run_batch_Fig5.m](fig5/run_batch_Fig5.m) reads each row of `JM_states` from `Fig5_constants.mat` and passes the corresponding fixed `J` and `M` to [run_NO_PSW.m](fig5/run_NO_PSW.m). Each run uses one of the three states `(0.5, 0.5)`, `(1.5, 0.5)`, and `(1.5, 1.5)` for every particle.

## MATLAB toolboxes

- **Parallel Computing Toolbox** is required. The simulations use `parpool`, `parfor`, and `spmd`. [Parallel pools](https://www.mathworks.com/help/parallel-computing/parpool.html) · [GPU computing](https://www.mathworks.com/help/parallel-computing/parallel.gpu.gpudevice.html)
- **Statistics and Machine Learning Toolbox** is required because particle initialization uses `randsample` to sample rotational states. The toolbox also provides `hist3` for optional image processing in step 8 and `Image_code_voigt.m`. [randsample](https://www.mathworks.com/help/stats/randsample.html) · [hist3](https://www.mathworks.com/help/stats/hist3.html)

## Tested environment and reproducibility

Reduced-size tests of the figure simulations were run on Windows with MATLAB R2026a.
The simulations sample initial particle positions, velocities, and laser fluctuations randomly. The scripts do not set a fixed random seed, so individual particle trajectories and velocities can differ between runs with the same settings.
