function run_batch_figS3(override)
% Run the figS3 PSW simulation for each quantum state used in the figure.

    arguments
        override (1, 1) struct = struct()
    end

    base_label = 'figS3';
    if isfield(override, 'label') && ~isempty(char(override.label))
        base_label = char(override.label);
    end

    states = {'00', '20', '21', '22', '40', '41', '42', '43', '44'};
    for state_idx = 1:numel(states)
        state = states{state_idx};
        override.jm_data = sprintf('alpha_JM_%s.txt', state);
        override.label = sprintf('%s_%s', base_label, state);
        fprintf('\nState %s (%d/%d)\n', state, state_idx, numel(states));
        run_all_figS3(override);
    end
end
