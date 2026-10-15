%% Energy Efficiency Calculation UPLINK

function [EnergyEfficiency] = EnergyEfficiencyCalc (SumCapacity, Performance, Power, PowerDiss)

EnergyEfficiency = (SumCapacity) / (sum(Performance .* Power) + sum(PowerDiss));

end