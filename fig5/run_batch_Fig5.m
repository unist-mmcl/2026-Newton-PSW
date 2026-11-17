% Run the Fig. 5 NO simulation for each fixed quantum state.

fig5_dir = fileparts(mfilename('fullpath'));
parameters = load(fullfile(fig5_dir, 'Fig5_constants.mat'));
states = parameters.JM_states;
validateattributes(states, {'numeric'}, {'2d', 'ncols', 2, 'nonempty', 'real', 'finite'});
parameters = rmfield(parameters, 'JM_states');

for state_idx = 1:size(states, 1)
    parameters.J = states(state_idx, 1);
    parameters.M = states(state_idx, 2);
    fprintf('\nFig. 5 state %d/%d (J=%.1f, M=%.1f)\n', ...
        state_idx, size(states, 1), parameters.J, parameters.M);
    run_NO_PSW(parameters);
end
