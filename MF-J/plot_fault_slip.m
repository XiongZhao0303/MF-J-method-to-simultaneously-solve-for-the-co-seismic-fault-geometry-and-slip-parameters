function plot_fault_slip(Fault, s, subdivisions,colormapInput)
    % plot_fault_slip Generates a 3D plot of fault slip distribution and 2D subplots for each fault.
    %
    % Input parameters:
    %   - Fault: A matrix containing the coordinate information of the fault surface vertices.
    %   - s: A column vector containing the components of the slip vector.
    %   - subdivisions: A cell array where each cell contains the number of subdivisions for each fault.
    %
    % Example:
    load mycolormap;

    % Determine the number of faults based on the length of subdivisions
    num_faults = length(subdivisions);
    
    % Initialize S1 and S2 as the two components of the slip vector
    S1 = s(1:2:end); % X-component of slip
    S2 = s(2:2:end); % Y-component of slip

    % Calculate slip size
    slip = sqrt(S1.^2 + S2.^2);

    % Extract coordinates of the fault surface vertices
    xfault = Fault(:, 4:7);
    yfault = Fault(:, 8:11);
    zfault = Fault(:, 12:15);

    % Create a 3D plot of the slip distribution on the fault surface
    figure;
    for i = 1:size(xfault, 1)
        patch(yfault(i, :), xfault(i, :), zfault(i, :), slip(i));
        hold on;
    end

    view(3); % Set 3D view
    grid on; % Show grid
    axis equal; % Set equal axis ratios
    rotate3d; % Enable rotation feature
    
    % Set colormap and colorbar
    colormap(flipud(colormapInput)); 
    cb = colorbar; 
    ylabel(cb, 'Slip / m'); % Set the label for the colorbar
    xlabel('Easting / km');
    ylabel('Northing / km');
    zlabel('Depth / km');
    set(gca, 'CameraPosition', [100 100 10]);
    
    % Create 2D subplots for each fault
    figure;
    startIndex = 1; % Initialize the starting index for patches
    
    for i = 1:num_faults
        % Calculate the total number of subdivisions for the current fault
        num_subdivisions = subdivisions{i}(1) * subdivisions{i}(2);
        
        % Create a new subplot for each fault
        subplot(num_faults, 1, i);
        
        % Loop through the subdivisions for the current fault
        for j = startIndex:(startIndex + num_subdivisions - 1)
            % Use the patch function to draw the fault surface in 2D
            patch(yfault(j, :), zfault(j, :), slip(j));
            hold on;
            
            % Plot arrows using quiver to show the slip vector components
            quiver(Fault(j, 2), Fault(j, 3), S1(j), S2(j), 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);
        end
        
        % Set colormap and colorbar for each subplot
        colormap(flipud(colormapInput)); 
        cb = colorbar; 
        ylabel(cb, 'Slip / m'); % Set the label for the colorbar
        caxis([min(slip), max(slip)]); % Set color axis limits

        % Title and labels for the current subplot
        title(['Fault ', num2str(i)]); 
        xlabel('Easting / km');
        zlabel('Depth / km');    
        
        % Update the starting index for the next fault
        startIndex = startIndex + num_subdivisions;
    end
end
