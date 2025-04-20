% File: simulate_load_shedding.m

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
