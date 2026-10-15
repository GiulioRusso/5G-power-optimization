%% Sum Capacity Calculation for UPLINK

function [SumCapacity, Capacity, SignalNoiseRatio, ChannelState] = SumCapacityCalc(h, PowerNoise, Power, Band, Interference)

[K,~] = size(Power);
SignalNoiseRatio = zeros(K, 1);
Capacity = zeros(K, 1);

% VEC. (K x 1) Channel State Information without Interference
ChannelState = zeros(K, 1);
for k=1:K
    % MAT. Linear Filter c (N x K) = h
    % MAT. c hermitian (K x N) = h hermitian
    % VEC. ck  = column k of h (N x 1) then hermitian
    ck = h(:,k)';
    % VEC. column k of h
    hk = h(:,k);
    % Num. of CSI
    Num = abs(ck * hk)^2;
    % Den. of CSI without Interference
    Den = PowerNoise * norm(ck)^2;
    
    % Interference terms (all elements of the column except the k one)
    if Interference == true
        for j=1:K
            if j~=k
                hj = h(:,j);
                Den = Den + Power(j) * abs(ck * hj)^2;
            end
        end
    end
    
    % CSI for user k
    ChannelState(k) = Num / Den;
    
    % VEC. (K x 1) SNR vector
    SignalNoiseRatio(k) = ChannelState(k) * Power(k);
    % VEC. (K x 1) Channel capacity for user number k
    Capacity(k) = Band * log2(1 + SignalNoiseRatio(k));
end

% SCA. Sum of all K Channel Capacities for the value i of Power
SumCapacity = sum(Capacity);

end