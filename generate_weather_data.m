% File: generate_weather_data.m

function weatherData = generate_weather_data()
    rng(0);  % Initialize random number generator for reproducibility

    % Generate synthetic solar irradiance data (W/m^2)
    hours = 0:23;
    solarIrradiance = zeros(1, 24);
    solarIrradiance(7:18) = max(0, 1000 * sin(pi * (hours(7:18) - 6) / 12));  % Sinusoidal pattern

    % Generate synthetic temperature data (°C)
    temperature = 15 + 10 * sin(pi * (hours - 6) / 12);

    % Save the weather data to a .mat file
    weatherData = struct('solarIrradiance', solarIrradiance, 'temperature', temperature);
    save('weatherData.mat', 'weatherData');
end
