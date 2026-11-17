function out_file = run_all_fig4(override)
% run_all_fig4(override) — PSW dispersion simulation for fig4
%
% Load fig4_constants.mat from this folder.
% Pass an optional struct to override individual constants.
%
% Run each sub-batch as a parfor task.
%   MAX=200,000 across 8 conditions yields 80 tasks of 20,000 particles each.
%   8 workers gave the best measured performance and saturated memory bandwidth.
%   The measured runtime was ~181 s (5.3× faster than the previous sequential run at ~960 s).

if nargin < 1 || isempty(override)
    override = struct();
end

tic

figure_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(figure_dir);
addpath(root_dir, fullfile(root_dir, 'sim'));

%% ── 0. Configure the parallel pool (set MaxNumWorkers) ───────────────────────
n_workers = 8;   % Best measured setting (8=181 s, 12=183 s, 16=187 s, 24=208 s)

% Restart the existing pool if its worker count differs from the requested count.
p = gcp('nocreate');
if ~isempty(p) && p.NumWorkers ~= n_workers
    delete(p);
    p = [];
end

if isempty(p)
    cluster = parcluster('local');
    if cluster.NumWorkers ~= n_workers
        cluster.NumWorkers = n_workers;  % Set MaxNumWorkers in the profile
        saveProfile(cluster);            % Save the profile changes
    end
    parpool(cluster, n_workers);
    fprintf('Started the parallel pool with %d workers\n', n_workers);
end

% Set 3 threads per worker to use 24 cores across 8 workers.
spmd
    maxNumCompThreads(3);
end


%% ── 1. Load constants ──────────────────────────────────────────────────────
batch_size = 5000;   % Default particle count per sub-batch
load(fullfile(figure_dir, 'fig4_constants.mat'));

% Apply the override struct fields to the local workspace.
% Handle label separately because it names the output file.
if isfield(override, 'label')
    result_label = override.label;
else
    result_label = '';
end
override_fields = fieldnames(override);
for fi = 1:numel(override_fields)
    if ~strcmp(override_fields{fi}, 'label')
        eval(sprintf('%s = override.(''%s'');', override_fields{fi}, override_fields{fi}));
    end
end
clear override_fields fi

%% ── 2. Image processing parameters ────────────────────────────────────────
rep_std       = 4.0;
rep_vel       = 2.7;
vel_pix       = 2.28;
atten_factorx = 1 / 0.9278;
atten_factory = 1 / 0.9453;

%% ── 3. Quantum state data file path ────────────────────────────────────────
jm_data_path = fullfile(root_dir, 'data', jm_data);
if ~isfile(jm_data_path)
    error('File not found: %s', jm_data_path);
end

%% ── 4. Preload linearly.mat (parfor broadcast variable) ─────────────────────
fprintf('[1/2] Loading linearly.mat...\n');
load(fullfile(root_dir, 'linearly.mat'));
[nj_lin, nm_lin, ~] = size(linearly);
linearly_2d = reshape(linearly, nj_lin*nm_lin, []);  % (nj*nm)×nI, shared with all workers
clear linearly;

%% ── 5. Precompute laser intensities ───────────────────────────────────────
max_I   = 2.0 / (pi * w0z * w0y) * sqrt(4.0 * log(2) / (pi * tau^2));
max_I   = max_I * atten1 * atten2 * 1e-3;

if numel(I1_list) ~= numel(I2_list)
    error('I1_list and I2_list must have the same number of elements (I1=%d, I2=%d)', numel(I1_list), numel(I2_list));
end

I1_vals = max_I * I1_list(:);   % n_groups×1
I2_vals = max_I * I2_list(:);   % n_groups×1

n_groups   = numel(I2_list);
n_sub      = ceil(MAX / batch_size);         % Sub-batches per condition (e.g., 10)
n_tasks    = n_groups * n_sub;               % Total parfor tasks (e.g., 80)

fprintf('[2/2] Starting time integration...\n');
fprintf('      Conditions: %d, sub-batches per condition: %d, particles per sub-batch: %d, total tasks: %d\n', ...
    n_groups, n_sub, batch_size, n_tasks);
fprintf('      With %d workers, expect ceil(%d/%d) = %d rounds\n', ...
    n_workers, n_tasks, n_workers, ceil(n_tasks/n_workers));

%% ── 6. Simulate sub-batches in parallel with parfor ────────────────────────
% Pack broadcast parameters into a struct for the helper function.
P.x_position_fwhm  = x_position_fwhm;
P.Dye_sigma         = Dye_sigma;
P.Dye_tau           = Dye_tau;
P.vx_caused_tilted  = vx_caused_tilted;
P.x_init_vel_sigma  = x_init_vel_sigma;
P.y_init_vel_sigma  = y_init_vel_sigma;
P.z_init_vel_off    = z_init_vel_off;
P.z_init_vel_sigma  = z_init_vel_sigma;
P.Ij_distribution   = Ij_distribution;
P.Ips_z             = Ips_z;
P.Ips_y             = Ips_y;
P.lambda            = lambda;
P.w0z               = w0z;
P.w0y               = w0y;
P.tau               = tau;
P.tstep             = tstep;
P.t_ini             = t_ini;
P.tend              = tend;

results = cell(n_tasks, 1);
particles_per_sub = min(batch_size, MAX - (0:n_sub-1) * batch_size);

parfor task = 1:n_tasks
    sub_idx = mod(task - 1, n_sub) + 1;
    results{task} = run_parfor_task(task, n_sub, particles_per_sub(sub_idx), jm_data_path, ...
        linearly_2d, nj_lin, nm_lin, I1_vals, I2_vals, P);
end

%% ── 7. Reassemble sub-batch results by condition ───────────────────────────
Cell_velocity = cell(n_groups, 1);
for g = 1:n_groups
    task_idx = (g-1)*n_sub + (1:n_sub);           % Task indices for this condition
    Cell_velocity{g} = vertcat(results{task_idx});
    results(task_idx) = {[]};
end
clear results

%% ── 8. Image processing (parfor) ───────────────────────────────────────────
% Uncomment the code below when image processing is needed. If the saved file is large, you can run this step separately later.
% fprintf('[3/3] Image processing (parfor)...\n');
% [img_matrix, sumx, sumy] = image_processing(Cell_velocity, I2_list, ...
%     rep_std, rep_vel, vel_pix, atten_factorx, atten_factory);

%% ── 9. Save results ────────────────────────────────────────────────────────
% Simulation controls
constants.n_workers     = n_workers;
constants.MAX           = MAX;
constants.batch_size    = batch_size;
% Image processing
constants.rep_std       = rep_std;
constants.rep_vel       = rep_vel;
constants.vel_pix       = vel_pix;
constants.atten_factorx = atten_factorx;
constants.atten_factory = atten_factory;
% Molecular beam initial conditions
constants.jm_data           = jm_data;
constants.x_position_fwhm   = x_position_fwhm;
constants.Dye_sigma          = Dye_sigma;
constants.Dye_tau            = Dye_tau;
constants.vx_caused_tilted   = vx_caused_tilted;
constants.x_init_vel_sigma   = x_init_vel_sigma;
constants.y_init_vel_sigma   = y_init_vel_sigma;
constants.z_init_vel_off     = z_init_vel_off;
constants.z_init_vel_sigma   = z_init_vel_sigma;
constants.Ij_distribution    = Ij_distribution;
constants.Ips_z              = Ips_z;
constants.Ips_y              = Ips_y;
% Laser initial conditions
constants.lambda = lambda;
constants.w0z    = w0z;
constants.w0y    = w0y;
constants.tau    = tau;
constants.atten1 = atten1;
constants.atten2 = atten2;
constants.tstep  = tstep;
constants.t_ini  = t_ini;
constants.tend   = tend;
constants.max_I = max_I;
if ~isempty(result_label)
    out_name = result_label;
else
    out_name = 'fig4';
end
result_dir = fullfile(figure_dir, 'result');
if ~isfolder(result_dir)
    mkdir(result_dir);
end
out_file = fullfile(result_dir, sprintf('result_%s.mat', out_name));
save(out_file, 'Cell_velocity', 'I1_list', 'I2_list', 'constants', '-v7.3');
fprintf('Results saved to %s\n', out_file);

fprintf('\n=== Simulation complete. '); toc

end
