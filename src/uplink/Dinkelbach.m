%% Dinkelbach UPLINK

function [lambda, PowerOptimal, SumCapacityOpt] = Dinkelbach (Band, ChannelState, Performance, PowerDiss, SumCapacity, PowerMax, h, PowerNoise, index, EnergyEfficiencyVec, lambda0, PowerOptimalMat)

[K,~] = size(ChannelState);

% SCA. Epsilon value for convergence
epsilon = 1e-3;

% SCA. Counter
n = 0;

% Energy Efficiency value. Once the peak is reached, we don't calculate any
% lambda after that
if index == 1
    lambda = lambda0;
else
    lambda = EnergyEfficiencyVec(index - 1, 1);
end

% VEC. Feasible Power values P_*
PowerFeasible = zeros(K, 1);

% VEC. Optimal Power values P_~
PowerOptimal = zeros(K, 1);

% SCA. Optimal Capacity value
SumCapacityOpt = SumCapacity;

% SCA. New objective function
F = (SumCapacity) - lambda * (sum(PowerMax .* Performance) + sum(PowerDiss));



while F > epsilon
    
    % Find xn_*
    for k=1:K
        PowerFeasible(k) = (Band / (lambda * Performance(k))) - (1 / ChannelState(k));
        % 0 < P_~ < PowerMax
        PowerOptimal(k) = max(min(PowerMax, PowerFeasible(k)), 0);
    end
    
    % Find f(xn_*)
    [SumCapacityOpt, C_EE, SINR_EE, ChannelState] = SumCapacityCalc(h, PowerNoise, PowerOptimal, Band, false);
    
    % F = f(xn_*) - lambda * g(xn_*)
    F = (SumCapacityOpt) - lambda * (sum(PowerOptimal .* Performance) + sum(PowerDiss));
    
    % lambda = f(xn_*) / g(xn_*)
    lambda = (SumCapacityOpt) / (sum(PowerOptimal .* Performance) + sum(PowerDiss));
    
    % n increase
    n = n + 1;
    
end



% if all the vector is 0 because we didn't entered the while loop
n = 0;
for k=1:K
    if PowerOptimal(k) == 0
        n = n + 1;
    end
end

% the optimal power allocation will be the previous
if index > 1 && n == K
    PowerOptimal = PowerOptimalMat(:,index - 1);
end



end