function velocity = time_integrate_NO(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend, molecule, use_gpu)
% TIME_INTEGRATE  Euler integration.

    if nargin < 18 || isempty(use_gpu)
        use_gpu = gpu_available();
    end

    if use_gpu && ~gpu_available()
        warning('NO_PSW:GpuUnavailable', 'GPU was requested, but no GPU is available. Falling back to CPU.');
        use_gpu = false;
    end

    if use_gpu
        g_dev = gpuDevice;
        nI = size(aligned_cost, 2);
        if size(aligned_cost, 1) == 1
            lookup_mode = 'shared lookup';
            fixed_lookup_bytes = nI * 8;
            bytes_per_particle = 1536;
        else
            lookup_mode = 'per-particle lookup';
            fixed_lookup_bytes = 0;
            bytes_per_particle = nI * 8 + 250;
        end
        memory_budget = max(1, g_dev.AvailableMemory * 0.8 - fixed_lookup_bytes);
        gpu_batch_size = floor(memory_budget / bytes_per_particle);
        gpu_batch_size = max(1, min(gpu_batch_size, MAX));

        n_batches = ceil(MAX / gpu_batch_size);
        fprintf('  GPU: %s | %.2f GB available | %s | %d particles/batch | %d batches\n', ...
            g_dev.Name, g_dev.AvailableMemory / 1e9, lookup_mode, gpu_batch_size, n_batches);

        gpu_tic = tic;
        for b = 1:n_batches
            idx = ((b - 1) * gpu_batch_size + 1):min(b * gpu_batch_size, MAX);
            b_sz = numel(idx);
            aligned_batch = slice_aligned_cost(aligned_cost, idx);
            velocity(idx, :) = gpu_batch_NO( ...
                position(idx, :), velocity(idx, :), Ij(idx), Ips(idx, :), ...
                I1_vec(idx), I2_vec(idx), aligned_batch, last_data(idx, :), ...
                lambda, w0z, w0y, tau, b_sz, tstep, t_ini, tend, molecule);
        end
        fprintf('  GPU integration completed in %.1f s\n', toc(gpu_tic));
    else
        t = t_ini;
        while t < tend
            position = position + velocity .* tstep;
            [~, last_data, a] = axay_opt_partial_NO( ...
                position, t + Ij, I1_vec, I2_vec, ...
                aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX, molecule);
            velocity = velocity + a .* tstep;
            t = t + tstep;
        end
    end
end

function velocity = gpu_batch_NO(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend, molecule)

    position = gpuArray(position);
    velocity = gpuArray(velocity);
    aligned_cost = gpuArray(aligned_cost);
    I1_vec = gpuArray(I1_vec);
    I2_vec = gpuArray(I2_vec);
    Ij = gpuArray(Ij);
    Ips = gpuArray(Ips);
    last_data = gpuArray(last_data);

    t = t_ini;
    while t < tend
        position = position + velocity .* tstep;
        [~, last_data, a] = axay_opt_partial_NO( ...
            position, t + Ij, I1_vec, I2_vec, ...
            aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX, molecule);
        velocity = velocity + a .* tstep;
        t = t + tstep;
    end
    velocity = gather(velocity);
end

function tf = gpu_available()
    try
        tf = gpuDeviceCount('available') > 0;
    catch
        tf = false;
    end
end

function aligned_slice = slice_aligned_cost(aligned_cost, idx)
    if size(aligned_cost, 1) == 1
        aligned_slice = aligned_cost;
    else
        aligned_slice = aligned_cost(idx, :);
    end
end

