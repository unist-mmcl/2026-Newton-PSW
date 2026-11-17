function P = NO_parameters()
%NO_PARAMETERS Default parameters for the NO PSW dispersion simulation.
%
% Edit J and M here, or pass overrides to run_NO_PSW:
%   run_NO_PSW(struct('J', 1.5, 'M', 0.5, 'MAX', 200000))

    %% Quantum state selected by the user
    P.J = 0.5;
    P.M = 0.5;

    %% Execution controls
    P.MAX = 20000;
    P.batch_size = 5000;
    P.use_parallel = true;
    P.n_workers = 8;
    P.threads_per_worker = 3;

 
 
    P.use_gpu = false;
    P.use_gpu_in_parfor = false;

    %% NO molecular constants
    P.molecule.name = 'NO';
    P.molecule.mass = 4.9826e-26;  % kg, 14N16O = 30.0061 amu
    P.molecule.polarizability.parallel = 2.529e-40;      % C m^2 V^-1
    P.molecule.polarizability.perpendicular = 1.60e-40;  % C m^2 V^-1

    %% Laser intensity conditions
    % A scalar I1_list is expanded to match vector I2_list.
    P.I1_list = 0.2;
    P.I2_list = [1.8; 3.2; 5.0; 7.2; 9.8; 12.8; 16.2; 20.0];

    %% Initial beam and molecular-beam parameters copied from the CS2 run
    P.x_position_fwhm = 0.000342;
    P.Dye_sigma = 3.55e-6;
    P.Dye_tau = 1.699e-9;
    P.vx_caused_tilted = -3.0;
    P.x_init_vel_sigma = 1.3418259;
    P.y_init_vel_sigma = 1.57112527;
    P.z_init_vel_off = 1732.0;
    P.z_init_vel_sigma = 10.04;
    P.Ij_distribution = 1000.0;
    P.Ips_z = 1642.0;
    P.Ips_y = 1822.0;

    %% IR pulse and beam parameters
    P.lambda = 1.064e-6;
    P.w0z = 2.0e-5;
    P.w0y = 2.0e-5;
    P.tau = 7.1e-9;
    P.atten1 = 0.959;
    P.atten2 = 0.948;
    P.tstep = 1.0e-11;
    P.t_ini = -3.0e-8;
    P.tend = 3.0e-8;

    %% Output
    P.result_label = '';
end
