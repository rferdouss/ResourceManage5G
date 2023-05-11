function [graph_link_capacity, U, h, M, NFV_Placed_in_XPU, used_edge_matrix_for_a_flow] = updateGraph(graph_link_capacity, shortestPath, flowdemand,servicechain, optimal_service_allocation_among_xpu, optimal_service_allocation_among_xpu_cores, h, M, U, NFV_Placed_in_XPU, used_edge_matrix_for_a_flow)
%% We need to update four things for the converted graph
% 1. link capacity needs to be decreased for original graph
% 2. VNF capacity needs to be decreased from original graph
% 3. VNF cost will be zero for next graph if VNF has already been placed
% 4. No other VNF can be placed on same core that is already being used

% max number of cores
%max_cores = size(M,3);

% no of services
%no_services = size(M,1);

% no of nodes
%V = size(M,2);

% demand of flow
df = flowdemand;%f{1}{3};

% service chain of flow
S_chain = servicechain;%f{1}{5}; 

%% 1. Graph Link capacity. Change the link capacities
y_cp = graph_link_capacity;
i=1;
j=1;
while (i<length(shortestPath))
    j=i+1;
    pathid1 = shortestPath(1,i);
    pathid2 = shortestPath(1,j);
    % subtract from the capacity of the edge in both directions
    y_cp(pathid1,pathid2) = y_cp(pathid1,pathid2)-df;
    y_cp(pathid2,pathid1) = y_cp(pathid2,pathid1)-df;
    
    %make the used edge matrix
    used_edge_matrix_for_a_flow(pathid1,pathid2) = used_edge_matrix_for_a_flow(pathid1,pathid2)+1;
    used_edge_matrix_for_a_flow(pathid2,pathid1) = used_edge_matrix_for_a_flow(pathid2,pathid1)+1;
    
    % increase the index value
    i =i+1;
end 
graph_link_capacity = y_cp;
%% 2, 3 and 4 Decrease VNF cost and no one else cnn use the same core to se
%{
for n = 1 : no_services
    for v = 1 : V
        for a = 1 : max_cores
            % service n is being used at core a of node v
         
            if XX(n,v,a) == 1
               % decrease the capacity of service on that node 
               U(n,v,a) = U(n,v,a) - df;  
               
               % no other service, other then current one, can be placed on same core
               M(:,v,a) = 0;
               M(n,v,a) = 1;
               
               % cost of putting service on same core/node is zero
               h(n,v,a) = 0;
               
            end
        end
    end
%}
for i =1: length(S_chain)
    service_id = S_chain(i,1);
    node_id = optimal_service_allocation_among_xpu(1,i);
    core_id = optimal_service_allocation_among_xpu_cores(1,i);
    
    % decrease the capacity of service on that node 
   % disp 'service bandwidth capacity in xpu and bandwidht demand of flow:';
   % service_id
   % node_id
    %core_id
    %U(service_id,node_id,core_id)
    %df
    U(service_id,node_id,core_id) = U(service_id,node_id,core_id) - df; 
    
    % no other service, can be placed on same core
    M(:,node_id,core_id) = 0;
    NFV_Placed_in_XPU(service_id,node_id,core_id) =1;
    
    % cost of putting service on same core/node is zero
    h(service_id,node_id,core_id) =0;
end 
end