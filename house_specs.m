% File: house_specs.m

function [voltage, mainCurrent, mainCircuitBreaker, allCircuits, allLoads, allLoadPriority] = house_specs()
    % House and Electrical Specifications
    voltage = 230;  % Volts
    mainCurrent = 100;  % Amperes
    mainCircuitBreaker = 100;  % Amperes

    % Define main and sub circuits
    allCircuits = {'LR1', 'LR2', 'DINING', 'BED1', 'BED2', 'BED3', 'BATH1', 'BATH2', 'GARAGE', 'KITCH1', 'KITCH2', 'LAUNDRY'};
    allLoads = [115, 138, 69, 92, 115, 115, 92, 92, 138, 207, 184, 161];  % Example loads in watts
    allLoadPriority = [1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0];  % Priority for load shedding
end
