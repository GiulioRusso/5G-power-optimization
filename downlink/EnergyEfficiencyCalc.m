%% Energy Efficiency Calculation DOWNLINK

function [EnergyEfficiency] = EnergyEfficiencyCalc (SumCapacity, Performance, PowerTotal, PowerDiss)

EnergyEfficiency = (SumCapacity) / (Performance .* PowerTotal + PowerDiss);

end