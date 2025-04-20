% File: solar_simulation.m

% Main script to run the solar energy simulation

% 1. Battery capacity (Wh)
batteryCapacity = input('Enter the battery capacity in Wh (e.g., 2000): ');

% 2. Number of solar panels
numSolarPanels = input('Enter the number of solar panels (e.g., 5): ');

% 3. Solar panel rating (W per panel)
solarPanelRating = input('Enter the solar panel rating in W (e.g., 300): ');

% 4. Load profiles input
disp('Please enter the load profiles as a matrix:');
disp('Rows represent different circuits, columns represent hours of the day (0 or 1)');
disp('Example input format for 12 circuits and 24 hours:');
disp('[1 1 0 0 ...; 0 1 0 0 ...; ...]');
loadProfiles = input('Enter the load profiles matrix (e.g., [1 1 0 0; 0 1 0 1; ...]): ');

% Set other values as constants
voltage = 230;  % Volts
mainCurrent = 100;  % Amperes
mainCircuitBreaker = 100;  % Amperes

% Define main circuits (you can change or expand this based on the house setup)
allCircuits = {'LR1', 'LR2', 'DINING', 'BED1', 'BED2', 'BED3', 'BATH1', 'BATH2', 'GARAGE', 'KITCH1', 'KITCH2', 'LAUNDRY'};

% Example load values for each circuit (W)
allLoads = [115, 138, 69, 92, 115, 115, 92, 92, 138, 207, 184, 161]; 

% Example priority values (1 for essential, 0 for non-essential)
allLoadPriority = [1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0]; 

% Calculate total solar power output from the user's input
totalSolarPower = solarPanelRating * numSolarPanels;

% Other fixed values related to the battery and solar system
batterySOC = 0.5;  % Initial state of charge (50%)
batteryEfficiency = 0.9;  % Battery efficiency (charge/discharge)
maxChargeRate = 1000;  % Max charge rate (W)
maxDischargeRate = 1000;  % Max discharge rate (W)

% Inverter efficiency curve as a function of load
inverterEfficiencyCurve = @(x) 0.95 - 0.05 * exp(-x/1000);

% Generate synthetic weather data (irradiance and temperature)
weatherData = generate_weather_data();

% Plot the generated weather data
plot_weather_data();

% Simulate load shedding with the provided parameters
[solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding(weatherData, voltage, allLoads, allLoadPriority, loadProfiles, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve);

% Visualize the results (shed loads, solar power, and battery state)
visualize_results(allCircuits, solarPowerOutput, batteryPower, batterySOCArray, shedLoads);


% --------- FUNCTION DEFINITIONS --------- %

function weatherData = generate_weather_data()
    rng(0);  % Initialize random number generator for reproducibility

    % Generate synthetic solar irradiance data (W/m^2)
    hours = 0:23;
    solarIrradiance = zeros(1, 24);
    solarIrradiance(7:18) = max(0, 1000 * sin(pi * (hours(7:18) - 6) / 12));  % Sinusoidal pattern

    % Generate synthetic temperature data (°C)
    temperature = 15 + 10 * sin(pi * (hours - 6) / 12);

    % Save the weather data to a structure
    weatherData = struct('solarIrradiance', solarIrradiance, 'temperature', temperature);
end

function plot_weather_data()
    % Load the generated weather data
    weatherData = generate_weather_data();

    % Plot the generated data for visualization
    figure;
    subplot(2, 1, 1);
    plot(0:23, weatherData.solarIrradiance);
    xlabel('Hour of Day');
    ylabel('Solar Irradiance (W/m^2)');
    title('Synthetic Solar Irradiance Data');
    grid on;

    subplot(2, 1, 2);
    plot(0:23, weatherData.temperature);
    xlabel('Hour of Day');
    ylabel('Temperature (°C)');
    title('Synthetic Temperature Data');
    grid on;
end

function [solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding(weatherData, voltage, allLoads, allLoadPriority, loadProfiles, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve)
    timeSteps = 24;  % Simulating for 24 hours
    solarIrradiance = weatherData.solarIrradiance;  % Solar irradiance data

    % Initialize variables for storing results
    solarPowerOutput = zeros(1, timeSteps);
    batteryPower = zeros(1, timeSteps);
    batterySOCArray = zeros(1, timeSteps);
    shedLoads = zeros(size(loadProfiles));

    % Simulation Loop
    for t = 1:timeSteps
        % Calculate solar power output based on irradiance
        solarPowerOutput(t) = totalSolarPower * (solarIrradiance(t) / 1000) * inverterEfficiencyCurve(solarIrradiance(t));

        % Initialize the result matrix
        dailyLoad = zeros(size(loadProfiles));

        % Multiply each element in allLoads with the corresponding row in loadProfiles
        for i = 1:length(allLoads)
            dailyLoad(i, :) = loadProfiles(i, :) * allLoads(i);
        end
        % Sum all columns to get total hourly loads
        totalHourlyLoads = sum(dailyLoad, 1);

        % Calculate the power balance
        powerBalance = solarPowerOutput(t) - totalHourlyLoads(t);

        % Battery charging and discharging logic
        if powerBalance > 0
            % Excess solar power available, charge the battery
            chargePower = min(powerBalance, maxChargeRate);
            chargeEnergy = chargePower * batteryEfficiency;
            batterySOC = min(batterySOC + chargeEnergy / batteryCapacity, 1);
            batteryPower(t) = chargePower;
        else
            % Insufficient solar power, discharge the battery
            dischargePower = min(abs(powerBalance), maxDischargeRate);
            dischargeEnergy = dischargePower / batteryEfficiency;
            if batterySOC * batteryCapacity >= dischargeEnergy
                batterySOC = batterySOC - dischargeEnergy / batteryCapacity;
                batteryPower(t) = -dischargePower;
            else
                % Battery capacity is not sufficient, perform load shedding
                deficit = abs(powerBalance) - batterySOC * batteryCapacity * batteryEfficiency;
                for i = 1:length(allLoads)
                    if allLoadPriority(i) == 0 && deficit > 0
                        shedLoads(i, t) = dailyLoad(i, t);
                        deficit = deficit - dailyLoad(i, t);
                        dailyLoad(i, t) = 0;
                    end
                end

                % If deficit still exists, shed priority loads
                if deficit > 0
                    for i = 1:length(allLoads)
                        if allLoadPriority(i) == 1 && deficit > 0
                            shedLoads(i, t) = dailyLoad(i, t);
                            deficit = deficit - dailyLoad(i, t);
                            dailyLoad(i, t) = 0;
                        end
                    end
                end

                batteryPower(t) = -batterySOC * batteryCapacity * batteryEfficiency;
                batterySOC = 0;
            end
        end

        % Store the state of charge
        batterySOCArray(t) = batterySOC;
    end
end

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
    grid on;

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
    grid on;

    % Add legend
    legend({'Solar Power Output'}, 'Location', 'northeastoutside');
    hold off;

    % Create a table to display simulation results
    create_results_table(solarPowerOutput, batteryPower, batterySOCArray, shedLoads, allCircuits);

end

function create_results_table(solarPowerOutput, batteryPower, batterySOCArray, shedLoads, allCircuits)

    % Create a time vector for the tablewer', batterySOCArray', shedLoadsSummary', ...
                         'VariableNames'; {'Hour', 'Solar_Power_Output_W', 'Battery_Power_W', 'Battery_SOC', 'Total_Shed_Loads_W'};

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
            'RowName', [], 'Position', [5 5 600 600]);

end
