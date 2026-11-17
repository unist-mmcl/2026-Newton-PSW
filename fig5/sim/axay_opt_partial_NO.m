function [cos2theta, last_data, acceleration, temp] = axay_opt_partial_NO( ...
        position, t, I1, I2, aligned_value, last_data, IR_position_stability, ...
        lambda, w0z, w0y, tau, MAX, molecule)
%AXAY_OPT_PARTIAL_NO Compute optical dipole acceleration for NO.

    epsilon0 = 8.854187e-12;
    c = 2.998e8;
    inv_mass = 1 / molecule.mass;
    alpha_parallel = molecule.polarizability.parallel;
    alpha_perpendicular = molecule.polarizability.perpendicular;

    gauss_t = exp(-4 * log(2) * t.^2 / tau^2);
    q = 2 * pi * position(:, 1) / lambda;

    I_fringe = I1 + I2 + 2 .* sqrt(I1 .* I2) .* cos(2 .* q);
    spatial_z = exp(-2 .* (position(:, 3) + IR_position_stability(:, 1)).^2 / w0z^2);
    spatial_y = exp(-2 .* (position(:, 2) + IR_position_stability(:, 2)).^2 / w0y^2);
    I = I_fringe .* spatial_z .* spatial_y;

    n_cols = size(aligned_value, 2);
    I_scaled = I .* gauss_t .* 1e-13;
    I_scaled = min(max(I_scaled, 0), n_cols - 1);

    I_floor = floor(I_scaled);
    frac = I_scaled - I_floor;
    upper_edge = I_floor >= n_cols - 1;
    I_floor(upper_edge) = n_cols - 2;
    frac(upper_edge) = 1;

    col_idx = I_floor + 1;
    if size(aligned_value, 1) == 1
        key = col_idx;
        lower_value = aligned_value(col_idx);
        upper_value = aligned_value(col_idx + 1);
        cos2theta = lower_value(:) + frac .* (upper_value(:) - lower_value(:));
    else
        row_idx = local_row_index(MAX, position);
        key = sub2ind(size(aligned_value), row_idx, col_idx);
        cos2theta = aligned_value(key) + frac .* (aligned_value(key + MAX) - aligned_value(key));
    end

    alpha = (alpha_parallel - alpha_perpendicular) .* cos2theta + alpha_perpendicular;
    mm = 0.5 * I / epsilon0 / c;

    valid_prev = last_data(:, 2) ~= 0;
    d_alpha_dI = (last_data(:, 1) - alpha) ./ (last_data(:, 2) - I .* gauss_t);
    d_alpha_dI = d_alpha_dI .* valid_prev;
    alpha_eff = alpha + I .* d_alpha_dI .* gauss_t;

    mask = (I ~= 0);
    acceleration = zeros(MAX, 3, 'like', position);
    n_k = 2 * pi / lambda;

    x_grad = 2 .* n_k .* (2 .* sqrt(I1 .* I2) .* sin(2 .* q)) ./ I_fringe;
    acceleration(mask, 1) = -mm(mask) .* x_grad(mask) .* inv_mass ...
        .* gauss_t(mask) .* alpha_eff(mask);

    acceleration(mask, 2) = -4 .* mm(mask) .* position(mask, 2) .* inv_mass ...
        .* gauss_t(mask) ./ w0y^2 .* alpha_eff(mask);

    acceleration(mask, 3) = -4 .* mm(mask) .* position(mask, 3) .* inv_mass ...
        .* gauss_t(mask) ./ w0z^2 .* alpha_eff(mask);

    last_data(:, 1) = alpha;
    last_data(:, 2) = I .* gauss_t;

    if nargout >= 4
        row_idx = local_row_index(MAX, position);
        temp = [I, I_scaled, I_floor, row_idx, col_idx, key, cos2theta];
    end
end

function row_idx = local_row_index(MAX, like_array)
    row_idx = (1:MAX)';
    if isa(like_array, 'gpuArray')
        row_idx = gpuArray(row_idx);
    end
end

