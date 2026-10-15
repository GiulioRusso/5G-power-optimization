%% Sum Capacity Calculation for DOWNLINK

function [SumCapacity, Capacity, SignalNoiseRatio, ChannelState] = SumCapacityCalc(h, PowerNoise, Power, Band, Interference)

[K,~] = size(Power);
SignalNoiseRatio = zeros(K, 1);
Capacity = zeros(K, 1);

% VEC. (K x 1) Channel State Information without Interference
ChannelState = zeros(K, 1);
for k=1:K
    % VEC. hk  = column k of h (N x 1)
    hk = h(:,k); 
    % VEC. column k of beamforming vector q
    qk = hk ./ norm(hk);
    % Num. of CSI
    Num = abs(hk' * qk)^2;
    % Den. of CSI without Interference
    Den = PowerNoise;
    
    % Interference terms (all elements of the column except the k one)
    if Interference == true
        for j=1:K
            if j~=k
                hj = h(:,j);
                qj = hj ./ norm(hj);
                Den = Den + Power(j) * abs(hk' * qj)^2;
            end
        end
    end
    
    % CSI for the user k
    ChannelState(k) = Num / Den;
    
    % VEC. (K x 1) SNR vector
    SignalNoiseRatio(k) = ChannelState(k) * Power(k);
    % VEC. (K x 1) Channel capacity for user number k
    Capacity(k) = Band * log2(1 + SignalNoiseRatio(k));
end

% SCA. Sum of all K Channel Capacities for the value i of Power
SumCapacity = sum(Capacity);

end