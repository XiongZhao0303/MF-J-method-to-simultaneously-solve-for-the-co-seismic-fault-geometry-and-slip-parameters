function plot_scatter_data(lose, def, d_o, lon,lat)
    % Plot scatter data based on the change of sign in lose
    % Inputs:
    %   losn - Array of LOS values
    %   def - Deformation values
    %   d_o - Observed deformation values
    %   XY - Coordinates matrix [Easting, Northing]

    % Identify the index where losn changes sign
    j = find_sign_change_index(lose);



    % Plotting
    figure;
    % First subplot
    subplot(2,3,1)
    scatter(lon(1:j), lat(1:j), 15, def(1:j),'filled');
    caxis([-max(def(1:j)), max(def(1:j))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on

    % Second subplot
    subplot(2,3,2)
    scatter(lon(1:j), lat(1:j), 15, d_o(1:j),'filled');
    caxis([-max(def(1:j)), max(def(1:j))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on

    % Third subplot
    subplot(2,3,3)
    scatter(lon(1:j), lat(1:j), 15, d_o(1:j) - def(1:j),'filled');
    caxis([-max(def(1:j)), max(def(1:j))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on

    % Fourth subplot
    subplot(2,3,4)
    scatter(lon(j+1:end), lat(j+1:end), 15, def(j+1:end),'filled');
    caxis([-max(def(j+1:end)), max(def(j+1:end))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on

    % Fifth subplot
    subplot(2,3,5)
    scatter(lon(j+1:end), lat(j+1:end), 15, d_o(j+1:end),'filled');
    caxis([-max(def(j+1:end)), max(def(j+1:end))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on

    % Sixth subplot
    subplot(2,3,6)
    scatter(lon(j+1:end), lat(j+1:end), 15, d_o(j+1:end) - def(j+1:end),'filled');
    caxis([-max(def(j+1:end)), max(def(j+1:end))]);
    colorbar;
    colormap(jet)
    xlabel('Easting /km'); ylabel('Northing /km');
    hold on
end

function j = find_sign_change_index(lose)
    % Find the index where losn changes sign
    for i = 1:length(lose)-1
        if lose(i) * lose(i+1) < 0
            j = i + 1;
            return;
        end
    end
    j = length(lose); % Return the last index if no sign change is found
end
