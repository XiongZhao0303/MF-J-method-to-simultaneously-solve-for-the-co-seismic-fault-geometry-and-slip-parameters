function plot_convergence_para(para_acc, fault_count,num_eachF_para)
    % plot_convergence_para Generates plots of convergence parameters for fault analysis.
    %
    % Input parameters:
    %   - para_acc: A matrix containing parameter values for convergence analysis.
    %   - fault_count: The number of faults (1 or 2).
    %
    % This function plots various parameters related to fault convergence, including
    % the components of the slip vector and the strike and dip angles for each fault.
    
    % Define legends for the plots
    legends = cell(1, fault_count*num_eachF_para+2); % Preallocate the legend array for up to 12 entries
    for i = 1:fault_count
        legends{(i-1)*num_eachF_para + 1} = sprintf('EX%d', i);     % X-component of slip for fault i
        legends{(i-1)*num_eachF_para + 2} = sprintf('EY%d', i);     % Y-component of slip for fault i
        legends{(i-1)*num_eachF_para + 3} = sprintf('EZ%d', i);     % Z-component of slip for fault i
        legends{(i-1)*num_eachF_para + 4} = sprintf('Strike%d', i); % Strike angle for fault i
        legends{(i-1)*num_eachF_para + 5} = sprintf('Dip%d', i);    % Dip angle for fault i
    end
    legends{num_eachF_para * fault_count + 1} = 'Sig_d';          % Variance component
    legends{num_eachF_para * fault_count + 2} = 'Sig_alpha';      % Variance component
 % Set up the subplot layout
    num_plots = fault_count * num_eachF_para + 2; % Total number of subplots
    rows = ceil(num_plots / 4);       % Calculate the number of rows
    cols = min(4, num_plots);         % Maximum of 4 columns
    % Set up the subplot layout for the plots
    figure;

    % Loop to plot each parameter graph
    for i = 1:(5 * fault_count+2) % Loop based on the number of faults
        subplot(rows, cols, i); % Create a subplot
        plot(para_acc(:, num_eachF_para * fault_count + 3), para_acc(:, i), 'c', 'LineWidth', 2); % Plot data
        legend(legends{i}); % Add legend for the current plot
        hold on; % Keep the current plot active for overlays
    end
end
