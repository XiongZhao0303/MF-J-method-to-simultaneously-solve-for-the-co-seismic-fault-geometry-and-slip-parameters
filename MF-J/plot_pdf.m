function plot_pdf(sample1, fault_count,num_eachF_para)
    % Plot PDF data for given samples.
    % Inputs:
    %   sample1 - Matrix containing sampling period data.
    %   fault_count - Number of faults 
    % Validate the input parameter for the number of faults
    % Define the legend labels for the plots
    legends = cell(1, fault_count*num_eachF_para+2); % Preallocate the legend array
    for i = 1:fault_count
        legends{(i-1)*num_eachF_para + 1} = sprintf('EX%d', i);     % Horizontal displacement for fault i
        legends{(i-1)*num_eachF_para + 2} = sprintf('EY%d', i);     % Vertical displacement for fault i
        legends{(i-1)*num_eachF_para + 3} = sprintf('EZ%d', i);     % Depth displacement for fault i
        legends{(i-1)*num_eachF_para + 4} = sprintf('Strike%d', i);  % Strike angle for fault i
        legends{(i-1)*num_eachF_para + 5} = sprintf('Dip%d', i);     % Dip angle for fault i
    end
    legends{fault_count * num_eachF_para + 1} = 'Sig_d';          % Variance component
    legends{fault_count * num_eachF_para + 2} = 'Sig_alpha';      % Variance component

    % Set up the subplot layout
    num_plots = fault_count * num_eachF_para + 2; % Total number of subplots
    rows = ceil(num_plots / 4);       % Calculate the number of rows
    cols = min(4, num_plots);         % Maximum of 4 columns
    figure;
    
    % Loop to plot each parameter histogram
    for i = 1:num_plots
        subplot(rows, cols, i);
        % Create a histogram with specified color properties suitable for scientific publications
        h = histogram(sample1(:, i), 'Normalization', 'count', ...
            'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'black', 'FaceAlpha', 0.7);
        
        % Obtain histogram values for further analysis (if needed)
        counts = h.Values;
        edges = h.BinEdges;
        width = 0.5 * h.BinWidth; % Calculate the width of the bins
        mid = edges(1:end-1) + width; % Calculate the midpoints of the bins
        [mm, nn] = max(counts); % Find the maximum count for additional analysis if required
        hold on;

        % Add a legend for the current plot
        legend(legends{i});
    end
end
