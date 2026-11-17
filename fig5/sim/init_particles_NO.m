function [particle, laser] = init_particles_NO(MAX, P)
%INIT_PARTICLES_NO Initialize NO particles for one selected J,M state.

    particle = zeros(MAX, 11);
    particle(:, 10) = P.J;
    particle(:, 11) = P.M;

    laser = zeros(MAX, 5);
    laser(:, 1) = P.Dye_sigma * randn(MAX, 1);
    laser(:, 2) = P.Dye_tau * randn(MAX, 1) / sqrt(3);
    laser(:, 3) = (-P.Ij_distribution + randi([0, round(P.Ij_distribution * 2)], MAX, 1)) * 1e-12;
    laser(:, 4) = P.Ips_z .* randn(MAX, 1) * 1e-9;
    laser(:, 5) = P.Ips_y .* randn(MAX, 1) * 1e-9;

    particle(:, 1) = P.x_position_fwhm / 2.355 * randn(MAX, 1);
    particle(:, 2) = P.Dye_sigma * randn(MAX, 1);

    particle(:, 4) = P.vx_caused_tilted + P.x_init_vel_sigma * randn(MAX, 1);
    particle(:, 5) = P.y_init_vel_sigma * randn(MAX, 1);
    particle(:, 6) = P.z_init_vel_off + P.z_init_vel_sigma * randn(MAX, 1);

    particle(:, 3) = (P.tend + laser(:, 2)) .* P.z_init_vel_off ...
        - (P.tend - P.t_ini + laser(:, 2)) .* particle(:, 6) ...
        + laser(:, 1);
end
