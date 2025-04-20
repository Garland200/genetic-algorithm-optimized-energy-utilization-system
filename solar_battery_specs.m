% File: solar_battery_specs.m

function [totalSolarPower, batteryCapacity, batterySOC, batteryEfficiency, maxChargeRate, maxDischargeRate, inverterEfficiencyCurve] = solar_battery_specs()
    % Solar System Specifications
    solarPanelRating = 300;  % Watts per panel
    numSolarPanels = 10;  % Number of solar panels
    totalSolarPower = solarPanelRating * numSolarPanels;  % Total power in watts

    % Battery Specifications
    batteryCapacity = 2000;  % Battery capacity in Wh
    batterySOC = 0.5;  % Initial state of charge (50%)
    batteryEfficiency = 0.9;  % Battery efficiency (both charge and discharge)
    maxChargeRate = 1000;  % Maximum charge rate in W
    maxDischargeRate = 1000;  % Maximum discharge rate in W

    % Inverter efficiency curve
    inverterEfficiencyCurve = @(x) 0.95 - 0.05 * exp(-x/1000);  % Example efficiency curve
end
