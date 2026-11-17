function velocity = time_integrate(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend, use_gpu)
% TIME_INTEGRATE  Euler integration.
%
% Optional use_gpu argument
%   true    — Use the GPU for standalone runs
%   false   — Use the CPU (recommended when called from a parfor worker)
%   omitted — Detect GPU availability automatically

    if nargin < 17
        use_gpu = gpuDeviceCount('available') > 0;
    end

    if use_gpu
        %% GPU path with automatic batch sizing based on available VRAM
        g_dev = gpuDevice;
        nI = size(aligned_cost, 2);
        bytes_per_particle = nI * 8 + 250;
        gpu_batch_size = floor(g_dev.AvailableMemory * 0.8 / bytes_per_particle);
        gpu_batch_size = min(gpu_batch_size, MAX);

        n_batches = ceil(MAX / gpu_batch_size);
        fprintf('  GPU: %s | %.2f GB available | up to %d particles per batch | %d batches\n', ...
            g_dev.Name, g_dev.AvailableMemory/1e9, gpu_batch_size, n_batches);

        tic;
        for b = 1:n_batches
            idx   = ((b-1)*gpu_batch_size + 1) : min(b*gpu_batch_size, MAX);
            b_sz  = numel(idx);
            velocity(idx,:) = gpu_batch( ...
                position(idx,:), velocity(idx,:), Ij(idx), Ips(idx,:), ...
                I1_vec(idx), I2_vec(idx), aligned_cost(idx,:), last_data(idx,:), ...
                lambda, w0z, w0y, tau, b_sz, tstep, t_ini, tend);
        end
        fprintf('  Completed in %.1f s\n', toc);

    else
        %% CPU path with multithreaded BLAS for parfor workers or systems without a GPU
        t = t_ini;
        while t < tend
            position = position + velocity .* tstep;
            [~, last_data, a, ~] = axay_opt_partial( ...
                position, t + Ij, I1_vec, I2_vec, ...
                aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX);
            velocity = velocity + a .* tstep;
            t = t + tstep;
        end
    end
end


function velocity = gpu_batch(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend)
    position     = gpuArray(position);
    velocity     = gpuArray(velocity);
    aligned_cost = gpuArray(aligned_cost);
    I1_vec       = gpuArray(I1_vec);
    I2_vec       = gpuArray(I2_vec);
    Ij           = gpuArray(Ij);
    Ips          = gpuArray(Ips);
    last_data    = gpuArray(last_data);

    t = t_ini;
    while t < tend
        position = position + velocity .* tstep;
        [~, last_data, a, ~] = axay_opt_partial( ...
            position, t + Ij, I1_vec, I2_vec, ...
            aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX);
        velocity = velocity + a .* tstep;
        t = t + tstep;
    end
    velocity = gather(velocity);
end
