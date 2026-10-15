%% Water Filling algorithm

function [PowerOpt, WaterLevel] = WaterFilling(ChannelState, PowerTotal)

% SCA. initial value of v that give us a sum Power grater than PowerTotal
v_i = 1e-5;
% SCA. bound of precision on PowerTotal
epsilon = 1e-5;

% Initialization of v and X
% v will be the value that give us a sum Power less than PowerTotal
v = v_i; 
% Power expressions from the solving of the Lagrangian Dual Function
X = max((1/v - 1./ChannelState), 0);
% Repeat this code until the sum of X is less than PowerTotal
while (sum(X) > PowerTotal)
    v = v * 1.001;
    X = max((1/v - 1./ChannelState), 0);
    % SumX = sum(X); % checking the sum of X during the process
end

% v is the value that give us a sum Power less than the total
% v_i is the value that give us a sum Power greater than the total

% v_av is the average between v_i and v
v_av = (v_i + v) / 2;
% Power values with v_av
X = max((1/v_av - 1./ChannelState), 0);

% v_i (sum(X) > PowerTotal) -------- v_av (sum(X) ~ PowerTotal) -------- v (sum(X) < PowerTotal)

% Repeat this code until the sum of X is into a range of PowerTotal
while ( (sum(X) < (PowerTotal - epsilon)) || (sum(X) > (PowerTotal + epsilon)) )
    
    % if v_av give us a sum of X less than the PowerTotal
    % we are between v_i and v_av
    if (sum(X) < (PowerTotal - epsilon))
        % update the upper value of v
        v = v_av;
    end
    
    % if v_av give us a sum of X greater than the PowerTotal
    % we are between v_av and v
    if (sum(X) > (PowerTotal + epsilon))
        % update the lower value of v
        v_i = v_av;
    end
    
    % update the v_av and the corrisponding Power values
    v_av = (v_i + v) / 2;
    X = max((1/v_av - 1./ChannelState), 0);
    
end

% The resulting v_av is the optimal v that give us the optimal Powers into
% a feasible range depending on which epsilon is set
PowerOpt = X;
% SumPowerOpt = sum(PowerOpt) % checking the sum of the PowerOpt
% Water Level of the Power values
WaterLevel = 1/v_av;

end