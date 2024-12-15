function [lb, ub] = slip_bound_constraints(f, M, rake, MAXSLIP)
    % ss: strike slip component, ds: dip slip component based on rake angle
    ss = abs(MAXSLIP * cosd(rake));  % Strike slip component
    ds = abs(MAXSLIP * sind(rake));  % Dip slip component

    % Initialize lower and upper bounds
    ssu = zeros(M / 2, 1) + ss;
    ssd = zeros(M / 2, 1) - ss;
    dsu = zeros(M / 2, 1) + ds;
    dsd = zeros(M / 2, 1) - 0.025;%If the earthquake is primarily thrust, then theoretically the slip (dsd) would be zero, so you only need to provide a very small value.
    lb = zeros(M, 1);
    ub = zeros(M, 1);
    
    lb(1:2:end) = ssd;
    lb(2:2:end) = dsd;
    ub(1:2:end) = ssu;
    ub(2:2:end) = dsu;

    % Initialize a list to collect all boundary constraint indices
    num = [];
    total_units = 0;  % This will track the total number of units processed so far

    % Loop through each fault in f to calculate boundary indices
    for faultIdx = 1:length(f)
        % Extract the number of strike and dip units for the current fault
        n_fault_strike = f{faultIdx}(1) - 1;  % Strike count for this fault
        n_fault_dip = f{faultIdx}(2) - 1;     % Dip count for this fault

        % Initialize boundary parts for this fault
        num_bottom = [];  % Bottom edge (down)
        num_top = [];     % Top edge (up)
        num_left = [];    % Left edge (left)
        num_right = [];   % Right edge (right)
        
        % Iterate over strike and dip units to calculate boundary indices
        for i = 1:n_fault_strike
            for j = 1:n_fault_dip
                % Bottom edge (when j == 1)
                if j == 1
                    num_bottom = [num_bottom; (j - 1) * n_fault_strike + i + total_units];
                end
                % Top edge (when j == n_fault_dip)
                if j == n_fault_dip
                    num_top = [num_top; (j - 1) * n_fault_strike + i + total_units];
                end
                % Left edge (when i == 1)
                if i == 1
                    num_left = [num_left; (j - 1) * n_fault_strike + i + total_units];
                end
                % Right edge (when i == n_fault_strike)
                if i == n_fault_strike
                    num_right = [num_right; (j - 1) * n_fault_strike + i + total_units];
                end
            end
        end

        % Combine boundary indices for this fault into the num array
        num = [num; num_bottom; num_top; num_left; num_right];

        % Update the total number of units processed so far
        total_units = total_units + n_fault_strike * n_fault_dip;
    end

    % Map each fault unit to the slip vector (strike and dip)
    num_slip = [];
    for i = 1:length(num)
        num_str = 2 * (num(i) - 1) + 1;  % Strike component index
        num_dip = 2 * (num(i) - 1) + 2;  % Dip component index
        num_slip0 = [num_str; num_dip];   % Add both strike and dip to slip constraints
        num_slip = [num_slip; num_slip0];
    end
    
    % Apply boundary conditions (set to 0 or a very small value)
    lb(num_slip) = 0;
    ub(num_slip) = 0.0001;
end
