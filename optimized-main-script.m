% Main script to run the solar energy simulation with Genetic Algorithm optimization

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

% Run Genetic Algorithm optimization for load shedding, passing all required data
optimize_load_shedding_with_ga(loadProfiles, allLoads, allLoadPriority, batteryCapacity, batterySOC, totalSolarPower, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve, weatherData, allCircuits);


% --------- FUNCTION DEFINITIONS --------- %

% Initialize the population with random load shedding schedules
function population = initialize_population(populationSize, numCircuits, numHours)
    population = cell(1, populationSize);
    for i = 1:populationSize
        % Each individual is a random matrix representing load shedding decisions
        population{i} = randi([0 1], numCircuits, numHours);  % 1 = shed, 0 = keep
    end
end

% Evaluate the fitness of each individual in the population
function fitnessScores = evaluate_fitness(population, loadProfiles, allLoads, allLoadPriority, batteryCapacity, batterySOC, totalSolarPower, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve, weatherData)
    numIndividuals = length(population);
    fitnessScores = zeros(1, numIndividuals);

    for i = 1:numIndividuals
        schedule = population{i};

        % Simulate load shedding with the individual's schedule
        [solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding_with_schedule(weatherData, allLoads, allLoadPriority, loadProfiles, schedule, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve);

        % Fitness function: Minimize total shed load, maximize power availability
        totalShedLoad = sum(shedLoads(:));  % Sum of all shed loads
        remainingPower = sum(solarPowerOutput + batteryPower);  % Total available power

        % Fitness score = (Remaining Power - Shed Load) + priority factor
        priorityFactor = sum(sum(schedule .* repmat(allLoadPriority', 1, 24)));  % Penalty for shedding critical loads
        fitnessScores(i) = remainingPower - totalShedLoad - priorityFactor * 10;  % Penalize more for shedding priority loads
    end
end

% Selection function to choose parents based on fitness scores
function selectedParents = selection(population, fitnessScores)
    totalFitness = sum(fitnessScores);
    probabilities = fitnessScores / totalFitness;  % Normalize fitness scores to probabilities

    selectedParents = cell(1, length(population));
    for i = 1:length(population)
        selectedParents{i} = population{roulette_wheel_selection(probabilities)};
    end
end

% Roulette wheel selection function
function selectedIndex = roulette_wheel_selection(probabilities)
    cumulativeProb = cumsum(probabilities);
    randomValue = rand;
    selectedIndex = find(cumulativeProb >= randomValue, 1);
end

% Crossover function to combine two parents into offspring
function offspring = crossover(parents, crossoverRate)
    numParents = length(parents);
    offspring = cell(1, numParents);

    for i = 1:2:numParents
        if rand < crossoverRate && i < numParents
            % Select two parents and perform crossover
            parent1 = parents{i};
            parent2 = parents{i+1};

            % Single-point crossover
            crossoverPoint = randi([1, size(parent1, 2)]);
            offspring{i} = [parent1(:, 1:crossoverPoint), parent2(:, crossoverPoint+1:end)];
            offspring{i+1} = [parent2(:, 1:crossoverPoint), parent1(:, crossoverPoint+1:end)];
        else
            % No crossover, pass parents directly
            offspring{i} = parents{i};
            if i < numParents
                offspring{i+1} = parents{i+1};
            end
        end
    end
end

% Mutation function to introduce random changes
function mutatedOffspring = mutation(offspring, mutationRate)
    numOffspring = length(offspring);
    mutatedOffspring = offspring;

    for i = 1:numOffspring
        for j = 1:numel(offspring{i})
            if rand < mutationRate
                mutatedOffspring{i}(j) = ~offspring{i}(j);  % Flip 0 to 1 or 1 to 0
            end
        end
    end
end

% Simulate load shedding with a specific schedule
function [solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding_with_schedule(weatherData, allLoads, allLoadPriority, loadProfiles, schedule, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve)
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

        % Apply load shedding schedule
        dailyLoad = loadProfiles .* (1 - schedule);  % 1 - schedule applies the shedding

        % Sum all columns to get total hourly loads, corrected for dimension mismatch
        totalHourlyLoads = sum(dailyLoad .* repmat(allLoads', 1, size(dailyLoad, 2)), 1);

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
                        shedLoads(i, t) = min(dailyLoad(i, t), deficit);  % Shed only up to the deficit amount
                        deficit = deficit - shedLoads(i, t);
                        dailyLoad(i, t) = dailyLoad(i, t) - shedLoads(i, t);  % Reduce the load by the amount shed
                    end
                end

                % If deficit still exists, shed priority loads
                if deficit > 0
                    for i = 1:length(allLoads)
                        if allLoadPriority(i) == 1 && deficit > 0
                            shedLoads(i, t) = min(dailyLoad(i, t), deficit);  % Shed only up to the deficit amount
                            deficit = deficit - shedLoads(i, t);
                            dailyLoad(i, t) = dailyLoad(i, t) - shedLoads(i, t);  % Reduce the load by the amount shed
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

% Generate synthetic weather data (irradiance and temperature)
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

% Plot weather data for visualization
function plot_weather_data()
    weatherData = generate_weather_data();

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

% Visualize the results of the simulation
% Visualize the results of the simulation
function visualize_results(allCircuits, solarPowerOutput, batteryPower, batterySOCArray, shedLoads)
    figure('Name', 'Solar System Optimization with Load Shedding');

    % Shed Loads Visualization (Magnitude of Load Shed in Watts)
    subplot(2, 1, 1);
    
    % Display the magnitude of shed loads in actual units (e.g., Watts)
    imagesc(shedLoads);  % Show the actual amount of load shed
    
    title('Load Shedding Over 24 Hours (Magnitude of Shed Load in Watts)');
    xlabel('Hour of Day');
    ylabel('Circuit Number');
    
    yticks(1:length(allCircuits));
    yticklabels(allCircuits);
    
    xticks(1:24);
    xticklabels(0:23);
    
    % Add color bar for the magnitude of shed loads
    colorbar;
    
    % Set colormap to reflect larger shed loads with darker shades
    colormap(flipud(gray));
    
    % Adjust the color axis to match the range of shed loads in actual watts
    caxis([0 max(shedLoads(:))]);  % Adjust color axis to the maximum load shed value
    
    % Ensure grid lines are visible for better clarity
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

    % Display additional results in a table
    create_results_table(solarPowerOutput, batteryPower, batterySOCArray, shedLoads, allCircuits);
end

% Create a table to display simulation results
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
            'RowName', [], 'Position', [5 5 600 600]);

end

% Genetic Algorithm for optimizing load shedding
function optimize_load_shedding_with_ga(loadProfiles, allLoads, allLoadPriority, batteryCapacity, batterySOC, totalSolarPower, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve, weatherData, allCircuits)
    % Parameters for the Genetic Algorithm
    populationSize = 50;  % Number of individuals in the population
    generations = 100;  % Number of generations to run the algorithm
    mutationRate = 0.1;  % Probability of mutation
    crossoverRate = 0.8;  % Probability of crossover

    % Initialize Population (Random load shedding schedules)
    population = initialize_population(populationSize, length(allLoads), 24);

    % Run Genetic Algorithm for 'generations' iterations
    for gen = 1:generations
        % Evaluate the fitness of each individual
        fitnessScores = evaluate_fitness(population, loadProfiles, allLoads, allLoadPriority, batteryCapacity, batterySOC, totalSolarPower, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve, weatherData);

        % Selection (Select parents based on fitness scores)
        selectedParents = selection(population, fitnessScores);

        % Crossover (Generate offspring through crossover of parents)
        offspring = crossover(selectedParents, crossoverRate);

        % Mutation (Mutate offspring with a given mutation rate)
        mutatedOffspring = mutation(offspring, mutationRate);

        % Form the new population
        population = mutatedOffspring;

        % Print the best fitness score in the current generation
        disp(['Generation ' num2str(gen) ': Best Fitness = ' num2str(max(fitnessScores))]);
    end

    % Get the best load shedding schedule after all generations
    [~, bestIndex] = max(fitnessScores);
    bestSchedule = population{bestIndex};

    % Simulate the final optimized load shedding schedule
    [solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding_with_schedule(weatherData, allLoads, allLoadPriority, loadProfiles, bestSchedule, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve);

    % Visualize the results
    visualize_results(allCircuits, solarPowerOutput, batteryPower, batterySOCArray, shedLoads);
end
