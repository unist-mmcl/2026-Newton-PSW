function out_file = run_NO_PSW(override)
%RUN_NO_PSW Run the NO partial standing-wave dispersion simulation.
%
% Examples:
%   run_NO_PSW()
%   run_NO_PSW(struct('J', 1.5, 'M', 0.5, 'MAX', 200000))
%   run_NO_PSW(struct('use_parallel', false, 'use_gpu', true, 'batch_size', 200000))

    if nargin < 1 || isempty(override)
        override = struct();
    end

    tic
    root_dir = fileparts(mfilename('fullpath'));
    addpath(root_dir, fullfile(root_dir, 'sim'));

    P = NO_parameters();
    P = apply_override(P, override);
    P = normalize_inputs(P);

    linearly_path = fullfile(root_dir, 'linearly_NO.mat');
    if ~isfile(linearly_path)
        error('NO alignment table not found: %s', linearly_path);
    end

    L = load(linearly_path, 'linearly', 'J_values', 'M_values', 'intensity_grid');
    J_values = L.J_values(:).';
    M_values = L.M_values(:).';
    P.J_index = find(abs(J_values - P.J) < 1e-9, 1);
    P.M_index = find(abs(M_values - P.M) < 1e-9, 1);

    if isempty(P.J_index)
        error('J=%.1f is not available. Available J values: %s', P.J, mat2str(J_values));
    end
    if isempty(P.M_index)
        error('M=%.1f is not available. Available M values: %s', P.M, mat2str(M_values));
    end
    if P.M > P.J
        error('Invalid NO state: |M|=%.1f must be <= J=%.1f.', P.M, P.J);
    end
    if isnan(L.linearly(P.J_index, P.M_index, 1))
        error('linearly_NO.mat has no alignment data for J=%.1f, M=%.1f.', P.J, P.M);
    end

    [nj_lin, nm_lin, ~] = size(L.linearly);
    P.linearly_row = sub2ind([nj_lin, nm_lin], P.J_index, P.M_index);
    linearly_2d = reshape(L.linearly, nj_lin * nm_lin, []);
    aligned_lookup = linearly_2d(P.linearly_row, :);
    clear L linearly_2d

    [I1_list, I2_list] = expand_intensity_lists(P.I1_list, P.I2_list);
    max_I = 2.0 / (pi * P.w0z * P.w0y) * sqrt(4.0 * log(2) / (pi * P.tau^2));
    max_I = max_I * P.atten1 * P.atten2 * 1e-3;

    I1_vals = max_I * I1_list(:);
    I2_vals = max_I * I2_list(:);

    n_groups = numel(I2_list);
    n_sub = ceil(P.MAX / P.batch_size);
    n_tasks = n_groups * n_sub;

    fprintf('NO PSW simulation: J=%.1f, M=%.1f, MAX=%d, batch_size=%d\n', ...
        P.J, P.M, P.MAX, P.batch_size);
    fprintf('Shared alignment lookup: 1 x %d intensity points\n', numel(aligned_lookup));
    fprintf('Conditions: %d, sub-batches per condition: %d, total tasks: %d\n', ...
        n_groups, n_sub, n_tasks);

    results = cell(n_tasks, 1);
    use_parallel_now = P.use_parallel && n_tasks > 1;
    if use_parallel_now
        use_parallel_now = start_parallel_pool(P.n_workers, P.threads_per_worker);
    end

    P.use_parallel = use_parallel_now;

    if use_parallel_now
        fprintf('Running with parfor (%d workers). GPU in parfor: %d\n', ...
            P.n_workers, P.use_gpu_in_parfor);
        parfor task = 1:n_tasks
            results{task} = run_NO_task(task, n_sub, P.batch_size, P.MAX, ...
                aligned_lookup, I1_vals, I2_vals, P);
        end
    else
        fprintf('Running serially. GPU: %d\n', P.use_gpu);
        for task = 1:n_tasks
            results{task} = run_NO_task(task, n_sub, P.batch_size, P.MAX, ...
                aligned_lookup, I1_vals, I2_vals, P);
        end
    end

    Cell_velocity = cell(n_groups, 1);
    for g = 1:n_groups
        task_idx = (g - 1) * n_sub + (1:n_sub);
        Cell_velocity{g} = vertcat(results{task_idx});
    end

    constants = P;
    constants.max_I = max_I;
    constants.linearly_file = linearly_path;
    constants.J_values = J_values;
    constants.M_values = M_values;
    constants.shared_alignment_lookup = true;
    constants.n_alignment_points = numel(aligned_lookup);
    constants.created_by = mfilename;

    label = result_label(P);
    result_dir = fullfile(root_dir, 'result');
    if ~isfolder(result_dir)
        mkdir(result_dir);
    end
    out_file = fullfile(result_dir, sprintf('result_%s.mat', label));
    save(out_file, 'Cell_velocity', 'I1_list', 'I2_list', 'constants', '-v7.3');
    fprintf('Saved: %s\n', out_file);
    fprintf('Done. ');
    toc
end

function P = apply_override(P, override)
    names = fieldnames(override);
    for idx = 1:numel(names)
        name = names{idx};
        if strcmp(name, 'label')
            P.result_label = override.(name);
        else
            P.(name) = override.(name);
        end
    end
end

function P = normalize_inputs(P)
    P.J = double(P.J);
    P.M = abs(double(P.M));
    P.MAX = round(double(P.MAX));
    P.batch_size = round(double(P.batch_size));
    P.n_workers = round(double(P.n_workers));
    P.threads_per_worker = round(double(P.threads_per_worker));

    if ~isscalar(P.J) || ~isscalar(P.M)
        error('J and M must be scalar values.');
    end
    if P.MAX <= 0 || P.batch_size <= 0
        error('MAX and batch_size must be positive.');
    end
    if P.n_workers <= 0
        error('n_workers must be positive.');
    end
end

function [I1_list, I2_list] = expand_intensity_lists(I1_list, I2_list)
    I1_list = I1_list(:);
    I2_list = I2_list(:);
    if isscalar(I1_list) && ~isscalar(I2_list)
        I1_list = repmat(I1_list, size(I2_list));
    elseif isscalar(I2_list) && ~isscalar(I1_list)
        I2_list = repmat(I2_list, size(I1_list));
    elseif numel(I1_list) ~= numel(I2_list)
        error('I1_list and I2_list must have the same length, unless one is scalar.');
    end
end

function ok = start_parallel_pool(n_workers, threads_per_worker)
    ok = true;
    try
        p = gcp('nocreate');
        if ~isempty(p) && p.NumWorkers ~= n_workers
            delete(p);
            p = [];
        end
        if isempty(p)
            parpool('local', n_workers);
        end
        if threads_per_worker > 0
            spmd
                maxNumCompThreads(threads_per_worker);
            end
        end
    catch ME
        warning('NO_PSW:ParallelUnavailable', ...
            'Parallel setup failed (%s). Falling back to serial.', ME.message);
        ok = false;
    end
end

function label = result_label(P)
    if ~isempty(P.result_label)
        label = char(P.result_label);
    else
        label = sprintf('NO_J%s_M%s', half_label(P.J), half_label(P.M));
    end
    label = regexprep(label, '[^A-Za-z0-9_.-]', '_');
end

function s = half_label(x)
    s = strrep(sprintf('%.1f', x), '.', 'p');
end

