function slip_results(Faults, num_faults, subdivisions, s_o)
    % Separate slip components and calculate total slip
    SS = s_o(1:2:end-1);  % Strike-slip component
    DS = s_o(2:2:end);    % Dip-slip component
    slip_total = sqrt(SS.^2 + DS.^2);  % Total slip magnitude

    Mw = zeros(num_faults, 1);  % Array to store the magnitude (Mw) of each fault
    Mo = zeros(num_faults, 1);  % Array to store the seismic moment (Mo) of each fault

    index_start = 1;  % Starting index for slip values of each fault

    % Calculate the seismic moment and magnitude for each fault
    for i = 1:num_faults
        % Determine the total number of subdivisions (grids) for the current fault
        num_subdivisions = subdivisions{i}(1) * subdivisions{i}(2);

        % Extract the slip values corresponding to the current fault
        slip_for_fault = slip_total(index_start:index_start + num_subdivisions - 1);
        
        % Extract fault parameters and convert them to numeric format
        fault_params = cell2mat(Faults);  % Convert to numeric matrix
        fault_params = fault_params (index_start:index_start + num_subdivisions - 1,:);
        area = fault_params(:, end);  % Extract relevant parameters (e.g., width or slip area)

        % Calculate seismic moment (Mo)
        Mo(i) = sum(26e9 .* slip_for_fault' .* area .* 1e6);

        % Calculate earthquake magnitude (Mw)
        Mw(i) = (2/3) * log10(Mo(i)) - 6.06;

        % Update starting index for the next fault
        index_start = index_start + num_subdivisions;

        % Display the results for the current fault
        fprintf('%s: %.4e\n', sprintf('Mo%d', i), Mo(i));  % Output Mo for the current fault
        fprintf('%s: %.4f\n', sprintf('Mw%d', i), Mw(i));  % Output Mw for the current fault
    end
end
