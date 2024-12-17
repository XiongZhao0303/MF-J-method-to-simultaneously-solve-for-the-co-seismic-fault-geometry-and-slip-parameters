function [Faults, subdivisions, M_total] = fault_geometry(num_faults, fault_lengths, fault_widths, x, fault_params,num_eachF_para)
    % faults, subdivisions, and total slip count
    % INPUT:
    %   num_faults   - Number of faults
    %   fault_lengths - Array of fault lengths
    %   fault_widths  - Array of fault widths
    %   x    -  (EX, EY, EZ, Strike, Dip)
    %   fault_params  - Cell array of subdivisions along strike and dip for each fault
    % OUTPUT:
    %   Faults        - Cell array of fault geometries
    %   subdivisions  - Cell array of subdivisions for each fault
    %   M_total       - Total number of slip elements

    Faults = cell(num_faults, 1);      % Cell array to store fault geometries
    subdivisions = cell(num_faults, 1); % Cell array for subdivisions
    M_total = 0;                       % Initialize total slip count

    for i = 1:num_faults
        % Generate fault geometry
        Faults{i} = simulatefault(fault_lengths(i), fault_widths(i), ...
            x((i-1)*num_eachF_para + 1), x((i-1)*num_eachF_para + 2), ...
            x((i-1)*num_eachF_para + 3), x((i-1)*num_eachF_para + 4), ...
            x((i-1)*num_eachF_para + 5), fault_params{i});

        % Calculate subdivisions for slip calculation
        n_strike = fault_params{i}(1) - 1; % Number of subdivisions along strike
        n_dip = fault_params{i}(2) - 1;   % Number of subdivisions along dip
        subdivisions{i} = [n_strike, n_dip];

        % Update total slip count
        M_total = M_total + 2 * n_strike * n_dip;
    end
end
