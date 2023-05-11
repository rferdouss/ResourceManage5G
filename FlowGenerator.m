function [ F ] = FlowGenerator( V, no_flows, size_M, demand, sd_demand, latency, sd_latency )
%% Generates the flows
% V = number of nodes
% no_flows = total number of flows generated
% size_M = total number of services in the network
% demand = mean capacity demand of the flow
% sd_demand = standard deviation of the capacity demand
% latency = mean latency demand of the flow
% st_latency = mean latency demand of the flow

%% Flow Generator
F_total = no_flows;
F={};
for f=1:F_total
    %get start node
    s = randi(V);
    % get end node t such that s ~= t 
    while true
        t = randi(V);
        if t ~= s
            break;
        end
    end
    while true
        d = normrnd(demand,sd_demand);
        l = normrnd(latency, sd_latency);
        if d > 0 && l > 0
            break;
        else
            disp('demand or latency has negative value, it should be positive value!');
        end
    end
    
    
    s_chain_length = 1+randi(size_M-1); % between 2-4 
    s_c = [];
    for i=1:s_chain_length
        while true
            ser = randi(size_M);
            if isempty(find(s_c == ser))
                s_c = [s_c; ser];
                break;
            end
        end
    end
    F{f}= {s,t,d,l,s_c};
end

end

