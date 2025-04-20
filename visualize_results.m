function visualize_results(allCircuits, solarPowerOutput, batteryPower, batterySOCArray, shedLoads)

    figure('Name', 'Solar System Optimization with Load Shedding');

    % Shed Loads Visualization
    subplot(2, 1, 1);
    imagesc(shedLoads);
    title('Load Shedding Over 24 Hours');
    xlabel('Hour of Day');
    ylabel('Circuit Number');
    yticks(1:length(allCircuits));
    yticklabels(allCircuits);
    xticks(1:24);
    xticklabels(0:23);
    colorbar;
    colormap(flipud(gray));
    
    grid on;  % Add grid for shed loads visualization

    % Solar Power Output Visualization
    subplot(2, 1, 2);
    hold on;
    title('Solar Power Output');
    xlabel('Time (Hours)');
    ylabel('Power (W)');
    xticks(1:24);
    xticklabels(0:23);
    
    % Plot solar power output
    plot(1:24, solarPowerOutput, 'g', 'DisplayName', 'Solar Power Output');
    
    % Add grid to solar power output
    grid on;  % Add grid for solar power output visualization

    % Add legend
    legend({'Solar Power Output'}, 'Location', 'northeastoutside');
    hold off;

    % Create a table to display simulation results
    create_results_table(solarPowerOutput, batteryPower, batterySOCArray, shedLoads, allCircuits);

end

function create_results_table(solarPowerOutput, batteryPower, batterySOCArray, shedLoads, allCircuits)

    % Create a time vector for the table
    timeHours = (0:23)';  % Column vector of hours

    % Create the shed loads summary for the table
    shedLoadsSummary = sum(shedLoads, 1);  % Sum of shed loads for each hour

    % Combine the data into a table
    resultsTable = table(timeHours, solarPowerOutput', batterySOCArray', shedLoadsSummary', ...
                         'VariableNames', {'Hour', 'Solar_Power_Output_W', 'Battery_SOC', 'Total_Shed_Loads_W'});

    % Create a new figure for the table
    figure('Name', 'Simulation Results Table');

    % Create the table in a uicontrol
    uitable('Data', resultsTable{:,:}, 'ColumnName', resultsTable.Properties.VariableNames, ...
            'RowName', [], 'Position', [20 20 600 300]);

end
