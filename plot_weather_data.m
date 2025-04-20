% File: plot_weather_data.m

function plot_weather_data()
    % Load the generated weather data
    load('weatherData.mat');

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
