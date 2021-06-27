%% Dinkelbach DOWNLINK

function [lambda, PowerOptimal, SumCapacityOpt] = Dinkelbach (Band, ChannelState, Performance, PowerDiss, SumCapacity, PowerMax, h, PowerNoise, index, EnergyEfficiencyVecInit, EnergyEfficiencyVec, PowerOptimalMat)

[K,~] = size(ChannelState);

% SCA. Counter
n = 0;

% SCA. Epsilon value for convergence
epsilon = 1e-10;

% VEC. Feasible Power values P_*
PowerFeasible = zeros(K, 1);

% VEC. Optimal Power values P_~
PowerOptimal = zeros(K, 1);

% Energy Efficiency value
lambda = EnergyEfficiencyVecInit(index, 1);
lambdaprev = 0;

% SCA. Optimal Capacity value
SumCapacityOpt = SumCapacity;

% SCA. New objective function
F = (SumCapacity) - lambda * (sum(PowerMax .* Performance) + sum(PowerDiss));



while lambda - lambdaprev > epsilon % alias F > epsilon
    
    % Find xn_*
    for k=1:K
        PowerFeasible(k) = (Band / (lambda * Performance)) - (1 / ChannelState(k));
        % P_~ > 0
        PowerOptimal(k) = max(PowerFeasible(k), 0);
    end
    
    % Find f(xn_*)
    [SumCapacityOpt, C_EE, SINR_EE, ChannelState] = SumCapacityCalc(h, PowerNoise, PowerOptimal, Band, false);
    
    % F = f(xn_*) - lambda * g(xn_*)
    F = (SumCapacityOpt) - lambda * (sum(PowerOptimal .* Performance) + sum(PowerDiss));
    
    % lambda of iteration n-1
    lambdaprev = lambda;
    
    % lambda = f(xn_*) / g(xn_*) of iteration n
    lambda = (SumCapacityOpt) / (sum(PowerOptimal .* Performance) + sum(PowerDiss));
    
    % n increase
    n = n + 1;
    
end

% Check for number issue: if the EE go down, I use the previous value
if index > 1 && lambda < EnergyEfficiencyVec(index - 1, 1)
    lambda = EnergyEfficiencyVec(index - 1);
    PowerOptimal = PowerOptimalMat(:,index - 1);
end



% Introducing of the sum constraint
if sum(PowerOptimal) > PowerMax
    
    % SCA. Counter
    n = 0;
    
    % Energy Efficiency value
    lambda = EnergyEfficiencyVecInit(index, 1);
    lambdaprev = 0;
    
    % SCA. New objective function
    F = (SumCapacity) - lambda * (sum(PowerMax .* Performance) + sum(PowerDiss));
    
    while lambda - lambdaprev > epsilon % alias F > epsilon
        
        % Power allocation
        [PowerOptimal, WaterLevel] = WaterFillingEE(ChannelState, PowerMax, Band, lambda, Performance);
        
        % Find f(xn_*)
        [SumCapacityOpt, C_EE, SINR_EE, ChannelState] = SumCapacityCalc(h, PowerNoise, PowerOptimal, Band, false);
        
        % F = f(xn_*) - lambda * g(xn_*)
        F = (SumCapacityOpt) - lambda * (sum(PowerOptimal .* Performance) + sum(PowerDiss));
        
        % lambda of iteration n-1
        lambdaprev = lambda;
        
        % lambda = f(xn_*) / g(xn_*) of iteration n
        lambda = (SumCapacityOpt) / (sum(PowerOptimal .* Performance) + sum(PowerDiss));
        
        % n increase
        n = n + 1;
    end
    
end



end