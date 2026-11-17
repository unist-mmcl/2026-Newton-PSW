function result = run_NO_task(task, n_sub, batch_size, max_particles, ...
        aligned_lookup, I1_vals, I2_vals, P)
%RUN_NO_TASK Execute one condition/sub-batch task.

    g_idx = ceil(task / n_sub);
    s_idx = mod(task - 1, n_sub) + 1;
    first_particle = (s_idx - 1) * batch_size + 1;
    this_batch_size = min(batch_size, max_particles - first_particle + 1);

    if this_batch_size <= 0
        result = zeros(0, 3);
        return
    end

    [particle_t, laser_t] = init_particles_NO(this_batch_size, P);

    aligned_cost_t = aligned_lookup;
    I1_t = repmat(I1_vals(g_idx), this_batch_size, 1);
    I2_t = repmat(I2_vals(g_idx), this_batch_size, 1);

    use_gpu = P.use_gpu;
    if P.use_parallel && ~P.use_gpu_in_parfor
        use_gpu = false;
    end

    result = time_integrate_NO( ...
        particle_t(:, 1:3), particle_t(:, 4:6), ...
        laser_t(:, 3), laser_t(:, 4:5), ...
        I1_t, I2_t, aligned_cost_t, zeros(this_batch_size, 2), ...
        P.lambda, P.w0z, P.w0y, P.tau, this_batch_size, ...
        P.tstep, P.t_ini, P.tend, P.molecule, use_gpu);
end

