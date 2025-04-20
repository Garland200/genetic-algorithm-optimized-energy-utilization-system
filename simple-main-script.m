

% Generate synthetic weather data
weatherData = generate_weather_data();

% Plot weather data
plot_weather_data();

% Load house and electrical specifications
[voltage, mainCurrent, mainCircuitBreaker, allCircuits, allLoads, allLoadPriority] = house_specs();

% Load solar and battery specifications
[totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve] = solar_battery_specs();

% Load dynamic load profiles
loadProfiles = load_profiles();

% Simulate load shedding
[solarPowerOutput, batteryPower, batterySOCArray, shedLoads] = simulate_load_shedding(weatherData, voltage, allLoads, allLoadPriority, loadProfiles, totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve);

% Visualize the results
visualize_results(allCircuits, solarPowerOutput, batteryPower, batterySOCArray, shedLoads);
