%%--------------------------------------------------------------------------------------------------------------------------------------------------
                                        %Resource Management Algorithm (RMA)
%-------------------------------------------------------------------------------------------------------------------------------------------------
%Heuristic solution to service (NFV) placement and flow allocation problem for perflow
%Date :  August 2016

% Special code for running on multiple graph

% Function Name : runnerFunction_RMA_Heuristic_MultiGraph( V_, number_of_simulation )
% Description  : Main function. This function allows running of multiple simulations.  
% INPUT Parameters : 
%==========================
%   (a)Number of nodes in the network, 
%   (b) frequency of the simulation running
% for example, runnerFunction_RMA_Heuristic_Perflow(10,50) % will collect 1 to 50 runs of the flow allocation for Node size 10. 

%OUTPUT :
%==========================

% INPUT FILE:  
% This script load an input file containing the following information required for flow allocation and service placement
%   (1) Network topology : connected graph with N number of nodes and E number of links. Links/edges can be of two technology (i. mm-wave link, and ii. ethernet link)
%   (2) Flow Set -  containing F flows. Each flow is defined with the following info : i. source node, ii. destination node, iii. bandwidth demand, iv. latency, and v. FlowGraph/service chain demand
%   (3) Differentiation of flows - Currently two types of flows are used i. video, and ii. audio. The input file contains, 
%           i. ratio of different flows present in the Flow Set.
%           ii. bandwidth and latency of each types of flow (e.g., mean demand, standard deviation of demand, mean latency)
%   (4) Information about the XPU nodes (Processing node) in the network,
%          i. service placement possibility, capacity, and cost in each XPU, (ii) core info in each XPU
%-------------------------------------------------------------------------------------------------------------------------------------------------

%%
function [  ] = runnerFunction_RMA_Heuristic_Standard_MultiGraph( V_, number_of_simulation)

number_of_graph =1;  % number of graph topologies for the simulation
no_flows = 100;       % number of flows  

%% GLOBAL MATRIX - to hold the final result (over n number of graphs each containing V_number of nodes in the graph topology, each graph simulation contains m run on f flows)

%percentage of link utilization - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Global_Used_Edge_Matrix= zeros(number_of_graph,20);             % PERCENGAGE of used link (mm-wave + ethernet)
Global_Used_MMWAVE_Edge_Matrix= zeros(number_of_graph,20);      % percentage of used mm-wave link 
Global_Used_Ethernet_Edge_Matrix= zeros(number_of_graph,20);    % percentage of used ethernet link

% average link capacity utilization statistics - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Final_Avg_Link_Capacity_Utilization_Matrix= zeros(number_of_graph,20);  % PERCENGAGE  of link (mm-wave + ethernet) utilization
Final_Avg_Link_Capacity_Utilization_MMWAVE_Matrix= zeros(number_of_graph,20);  % percentage of mm-wave link utilization
Final_Avg_Link_Capacity_Utilization_Ethernet_Matrix= zeros(number_of_graph,20);  % percentage of ethernet link utilization

% average nfv capacity utilization - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Final_Avg_NFV_Capacity_Utilization_Matrix = zeros(number_of_graph,20);

% statistics of flow allocation - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Unsuccessful_flow_statistics = zeros(number_of_graph,20);  % number of unsuccessfull flow

% statistics about average path length - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Average_Path_Length_Global = zeros(number_of_graph,20);  % average path length

% Delay Matrix - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Delay_Statistics = zeros(number_of_graph,20); 

% Cost matrix - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Cost_Matrix = zeros(number_of_graph,20);  % cost

% service/NFV chain request and fulfill statistics - global variable (over n topology with v number of nodes, each having m run on allocating f flows)
Service_Placement_Statistics = zeros(number_of_graph,20);   % Number of total placed services
Service_Reuse_Statistics = zeros(number_of_graph,20);       % Number of total reused services
Service_Demand_Statistics = zeros(number_of_graph,20);      % Number of total service demand

%% -----------------------------------------------------------------------------------------------------------------------------------------------
%               STARTING OF THE SIMULATION FOR THE MULTIPLE GRAPH
%-------------------------------------------------------------------------------------------------------------------------------------------------
for str =1 : number_of_graph
    disp '----------------------------------------------------------------------------';
    fprintf('Random Graph ID : %d\n',str);
    disp '-----------------------------------------------------------------------------';

%% declaration of input variable
%disp('Number of nodes:')
%disp(V_)


%% LOAD INPUT FILE (with pre-generated flow)
V = V_;             % Number of nodes e.g., [10 15 20 25 30]
filename ='';
filename = strcat('Input/Input_N_',num2str(V),'_F_',num2str(no_flows),'_i_',num2str(str),'.mat'); % input file with the generated flow, network toptology, capacity, latency, service allocation information
load(filename);   % load the input file
%fprintf('Input file name : %s \n', filename);

%% Generate Network Flows
% ratio of flows that are video flows
videoFlows = 0.70;

% flows for video and voice
no_flows_video = round(no_flows*0.7);
no_flows_voice = no_flows-no_flows_video;

% number of services in the network
M_size = 4;
% Load M where services can be places! Given the total number of XPU and XFE nodes. Load O, which is number of cores available on each node Load metrics ServiceCapacity, capacity that can be satisfied by a service number of services
%M_size = size(M,1);

% % FOR VIDEO TRAFFIC
% mean demand
d = 5; % MB/s
% standard deviation of demand
d_sd = 1;
% mean latency of flow, in seconds

l = 100/1000; %100ms
% sd of flow
l_sd = 20/1000; %20ms

% generate flows
F_video = FlowGenerator(V,no_flows_video,M_size, d, d_sd, l, l_sd);

% % FOR VOICE TRAFFIC
% mean demand
d = 1; % MB/s
% standard deviation of demand
d_sd = 0.2;
% mean latency of flow, in seconds
l = 10/1000; %10ms
% sd of flow
l_sd = 2/1000; %2ms

% generate flows
F_voice = FlowGenerator(V,no_flows_voice,M_size, d, d_sd, l, l_sd);

F = [F_video, F_voice];
% random permutation of flows- Shuffle flows
F = {F{randperm(length(F))}};

%% Make M where services can be places! Given the total number of XPU and XFE nodes
%% Metric O, which is number of cores available on each node
V_expanded = size(A_con{1,1},1);

ratio_XPU = 0.3;
% no of XPU nodes
no_XPU = floor(V*ratio_XPU);
% no of XFE nodes
no_XFE = V - no_XPU;
% no of cores available on XPU nodes
no_core_XPU = 4;

M = zeros(M_size,V_expanded, no_core_XPU);

O = zeros(1,V_expanded);

for i=1:no_XPU
    while true
        % select a node
        r = randi(V);
        % if node is not already selected as XPU node, select that node
        if sum(M(:,r))==0
            break;
        end
    end
    % node r can host all services
    M(:,r,:) = 1;
    O(1,r) = no_core_XPU;
end

%% Metrics U, load that can be satisfied by a service
% no of flows that can be served by each service
no_flow_eachVNF = 10000; % in Mb/s

U = zeros(M_size,V_expanded, no_core_XPU);
for i = 1:M_size
    U(i,:,:) = no_flow_eachVNF;
end

%% Create h(n,i), fixed cost of placing VNF of type n on node i
cost = 1000;
h = zeros(M_size, V_expanded, no_core_XPU) + cost;
%% Semi-Global Variables (that remain active till the end of all simulation run for a single topology/graph), e.g., 50 runs for 100 flows with n nodes  
k =30;   % NUMBER OF shortest path to investigate from source node to destination

%percentage of used edge statistics (after m simulation run for a single topology)
num_used_edge_after_all_simulation= zeros(1,20);
percentage_of_used_mmwave_edge = zeros(1,20);
percentage_of_used_ethernet_edge = zeros(1,20);

% link capacity utilization (after m simulation run for a single topology)
link_capacity_utilization = zeros(1,20);
mmwave_link_capacity_utilization = zeros(1,20);
ethernet_link_capacity_utilization = zeros(1,20);
% link_utilization_matrix = zeros(20, number_of_simulation);
% mmwave_link_utilization_matrix = zeros(20, number_of_simulation);
% eth_link_utilization_matrix = zeros(20, number_of_simulation);

% average NFV capacity utilization (after m simulation run for a single topology)
nfv_capacity_utilization= zeros(1,20);

%service statistics (after m simulation run for a single topology)
service_reused = zeros(1,20);
service_placed = zeros(1,20);
service_demand = zeros(1,20);

total_service_demand=0;

total_mmwave_link=0;
total_eth_link=0;

% flow allocation statistics (after m simulation run for a single topology)
unsuccessfull_flow = zeros(1,20);
% average path length (after m simulation run for a single topology)
path_length_avg = zeros(1,20);   

%average flow allocation delay (after m simulation run for a single topology)
avg_flow_allocation_delay = zeros(1,20);
%average cost (after m simulation run for a single topology)
avg_cost =zeros(1,20);

tic;
%% RUNNING SIMULATION FOR M TIMES ON A GRAPH TOPOLOGY (each time allocate 100 flows on V nodes)
for run = 1:number_of_simulation
 
 %fprintf('Simulation run ID : %d \n', run);
%% Input Graph from loaded input file
y = A_con{1,1}; % graph connectivity
total_number_of_edge = sum(sum(y));  % total number of edges

y_c =A_con{1,2};    % fixed cost of using link
y_d =A_con{1,3};    % dynamic cost of using link
y_cp =A_con{1,4};   % graph link capacity
graph_capacity = A_con{1,4}; % graph link capacity
max_link_capacity = A_con{1,4};  % graph maximum link capacity

y_l =A_con{1,5};    % link latency
N = size(y,1); % Number of nodes
%disp 'Number of nodes after graph expansion';
%N
%% INPUT Flow parameters
ServiceMatrix =  M;
ServiceCapacity = U;
ServiceCostMatrix = h;


%% LOCAL VARIABLES active for a single simulation run only (for all flow allocation for once)
NFV_Placed_in_XPU = zeros(M_size,V, no_core_XPU);  % keep track of all the placed NFVs/Services(in one simulation)
used_edge_matrix=zeros(V,V); % matrix to store the used edge information (in one simulation)

% statistics about NFV placement and reuse
number_of_total_placed_services = 0;   % number of total service placed (new services)(in one simulation)
number_of_total_reused_services = 0;   % number of total service reused(in one simulation)
total_requested_services =0;           % number of total service requested  (in one simulation)

path_length=0;      %path length(in one simulation)
parametric_cost=0;  %parametric cost (in one simulation)
number_of_unsuccessful_flow_allocation=0; %flow allocation statistics (in one simulation)
avg_delay = 0;  % average delay (in one simulation)

%% Initial preparation before running the loop of flow allocation 
%first prepare the cost matrix (in order to run the K-Shortest path algorithm)
% Implementation of the standard K-shortest path with Hop count
y_temp= y;  % giving the adjacency matrix
y_temp( y_temp==0 )=Inf;  % for the requirement of the k-shortest path calculation code all the cost that is zero are replaced by 'INF', indicating the absence of coxt/edge
% Randomization of the flow by shuffle all flows
F = {F{randperm(length(F))}};
   
%% RUNNING MAIN LOOP - ALLOCATION OF ALL FLOW AND SERVICE PLACEMENT -- FOR ONE SIMULATION RUN
 tic %start counting the time
 for i= 1: no_flows % 1st for loop
    % take one flow
    f = F(i);   
    % other info about the flow
    S = f{1}{1};          % start node
    D = f{1}{2};          % end node
    d_f = f{1}{3};        % flow demand
    l_f = f{1}{4};        % latency demand
    S_chain = f{1}{5};    % service chaining
    
 %   disp '========================================================================';
 %COMMENT   fprintf('FLOW NO : %d,  START NODE %d, END NODE %d \n', i, S, D);
 %	fprintf('\n');   
   
    %% Defining the variables  (local). active for a single flow
    optimal_shortestPaths=zeros(1,1);  % storing the path
    optimal_Costs=0;            % cost of flow allocation along the optimal path
    optimal_service_allocation_among_xpu =zeros (1,length(S_chain)); % service placement in the XPU nodes along the shortest path
    optimal_service_allocation_among_xpu_cores =zeros (1,length(S_chain)); % service placement in the XPU cores along the shortest path
    optimal_number_of_service_reuse=0;
    optimal_number_of_new_service_placed = 0;    
    flag_Optimalpath=0;
	used_edge_matrix_for_a_flow = zeros(V,V); % for keeping the used edge/link for this flow allocation
    
    %find the k shortest paths from source S to destination D , calling function 'kShortestPath'
    [shortestPaths, totalCosts] = kShortestPath(y_temp, S, D, k);
  
	%%    
    %Flow allocation - running over the k shortest path to find the optimum one with placing required/desired services
    if isempty(shortestPaths)
      %  fprintf('No shortest path available between these nodes\n\n');
    else
        for j = 1: length(shortestPaths)  % 2nd for loop (running over all the k shortest path to find the optimum one for allocation this flow)
            %Defining local variables
            service_allocation_among_xpu = zeros (1,length(S_chain)); % keep the id of the xpu from the array 'XPU_List' that contains the list of XPUs
            service_allocation_among_xpu_core = zeros (1,length(S_chain));
            number_of_service_reuse=0;
            flag_service_placement_possibility=0; % this flag indicates if the services can be placed along this shortest path or not
                     
            %% FILTER LEVEL 1 : CHECK GOODNESS OF THE SHORTEST PATH (Constrain 1: exists at least one XPU, and Constrain 2. availability of enough link bandwidth and low latency)
            [Flag_Path_Goodness, XPU_List] = Check_Path_Goodness(shortestPaths{j},O,graph_capacity, y_l, d_f, l_f);
        
            %% FILTER LEVEL 2: CHECK THE POSSIBILITY OF FLOW ALLOCATION IN THE PATH - if this shortest path is good (contains at least one XPU, enough link bandwidth, and low latency) then go for next level to check if this shortest path has the potential to fulfil the flow-graph
            if(Flag_Path_Goodness==1)
               % disp 'Filter 1 : Path is good. It contains at least one XPU. Have sufficient bandwidth and low latency- going to check more for the possibility of satisfying the flow-graph ';   %[flag_Optimalpath, path_cost, service_allocation_xpu, service_allocation_xpu_core] = Check_Flow_Allocation_Possibility(shortestPaths{j}, totalCosts(j), XPU_List, S_chain, no_core_XPU, NFV_Placed_in_XPU, U, d_f, ServiceCostMatrix, ServiceMatrix); %check the possibility of service placement and flow allocation in this shortest path         
                if (all(~NFV_Placed_in_XPU))   % if this is the first flow/ no service is placed yet/ no existing service  
                   % disp 'No NFV is placed yet / NFV_in_XPU matrix is empty';    
                    [flag_flowallocation, totalcost,service_allocation_xpu, service_allocation_xpu_core]= Check_ServicePlacement(totalCosts(j), XPU_List, S_chain, no_core_XPU, ServiceCapacity,d_f, ServiceCostMatrix, ServiceMatrix, service_allocation_among_xpu, service_allocation_among_xpu_core);
                    service_allocation_among_xpu = service_allocation_xpu;
                    service_allocation_among_xpu_core = service_allocation_xpu_core;                    
                    flag_service_placement_possibility = flag_flowallocation;                   
                else % next flow/ NFV_in_XPU matrix is not empty / go for checking if existing services can be reused   
                    % Check 1 : Check for best option - No new service placement - Use of peviously placed services   
                    [flag_service_placement_possibility, totalcost, number_of_service_reuse, service_allocation_xpu, service_allocation_xpu_core] = check_NFV_reuse_possibility(totalCosts(j),XPU_List, S_chain, no_core_XPU, NFV_Placed_in_XPU, ServiceCapacity, d_f,ServiceCostMatrix, ServiceMatrix); %check the possibility of service placement and flow allocation in this shortest path         
                end  %end of if (all(~NFV_in_XPU)) 
                
                %% Now check the result of the service placement possibility -  if this path is a candidate for optimal path
                 if(flag_service_placement_possibility==1) % service placement is possible in this path
                    if(optimal_Costs==0)
                        optimal_Costs = totalcost;
                        optimal_shortestPaths=zeros(1,length(shortestPaths{j}));
                        optimal_shortestPaths = shortestPaths{j};
                        optimal_service_allocation_among_xpu = service_allocation_xpu;
                        optimal_service_allocation_among_xpu_cores = service_allocation_xpu_core;
                        optimal_number_of_service_reuse = number_of_service_reuse;
                        optimal_number_of_new_service_placed = length(S_chain) - optimal_number_of_service_reuse; 
                    end
                    if((optimal_Costs>0) && (totalcost <= optimal_Costs))
                        % check the path with minimum residual capacity (path_capacity - flow_capacity_demand) and choose that as optimum
                        % if path_cost is the shortest then update
                        optimal_Costs = totalcost;
                        optimal_shortestPaths=zeros(1,length(shortestPaths{j}));
                        optimal_shortestPaths = shortestPaths{j};
                        optimal_service_allocation_among_xpu = service_allocation_xpu;
                        optimal_service_allocation_among_xpu_cores = service_allocation_xpu_core;
                        optimal_number_of_service_reuse = number_of_service_reuse;
                        optimal_number_of_new_service_placed = length(S_chain) - optimal_number_of_service_reuse;                        
                    end                    
                    flag_Optimalpath=1;
                 else %  if(flag_service_placement_possibility==0)
                    % COMMENT: disp 'Filter 1 : Service placement error - All Service placement is not possible in this path';
                 end 
      
            else
            end 
        end % end of 2nd for loop (traversing the set of shortest paths found for a flow from source to destination node)				
       
		% if any shortes path is found, for this 	
		if(flag_Optimalpath == 1)
            
              % update the capacity graph (A_con, F_i, ServiceMatrix, O, ServiceCapacity, ServiceCostMatrix, V_org)
              [graph_capacity, ServiceCapacity,ServiceCostMatrix,ServiceMatrix, NFV_Placed_in_XPU, used_edge_matrix_for_a_flow] = updateGraph(graph_capacity, optimal_shortestPaths, d_f, S_chain, optimal_service_allocation_among_xpu, optimal_service_allocation_among_xpu_cores, ServiceCostMatrix, ServiceMatrix,ServiceCapacity, NFV_Placed_in_XPU,used_edge_matrix_for_a_flow);  
                   % [A_con, ServiceCapacity, ServiceCostMatrix, ServiceMatrix] = updateGraph(A_con, F_i, link_used_matrix, service_placement_matrix, ServiceCapacity, ServiceCostMatrix, ServiceMatrix);  
                  % disp 'Filter 2 : Flow is served along the optimal path after checking all the K shortest path'; 
               %COMMENT    disp 'Flow is served'; 
                  
                  % update the final used_edge_link matrix
                  used_edge_matrix = used_edge_matrix + used_edge_matrix_for_a_flow;
                  number_of_total_placed_services = number_of_total_placed_services + optimal_number_of_new_service_placed ;
                  number_of_total_reused_services = number_of_total_reused_services + optimal_number_of_service_reuse;
                  total_requested_services =  total_requested_services + length(f{1}{5}); % sum of requested services   
                  
                  path_length = path_length+length(optimal_shortestPaths);
                  parametric_cost =parametric_cost+optimal_Costs; 
                  
                  d = used_edge_matrix.* y;
                  %delay = [delay sum(sum(d))];   
                  avg_delay =  avg_delay +  sum(sum(d));
                  %used_edge_matrix
        else
          %COMMENT  disp 'Filter 2 : This flow allocation is not possible - no optimal path found among the list of K shortest paths';
          %COMMENT  disp 'This flow allocation is not possible'; 
          disp '---------';
          fprintf('flow id = %d, simulation id =%d, graph id = %d\n', i, run, str);
          fprintf('source %d, destination %d, total shortest path %d \n',S,D, length(shortestPaths));
          number_of_unsuccessful_flow_allocation=number_of_unsuccessful_flow_allocation+1;            
		end % end of if(flag_successful_flow_allocation == 1)
    end % end of if shortest path list is not empty 
	
    
    %% save the result after each 5 flows
            if mod(i,5) == 0
            % calculate the ratio of used edges of different technology
            %% calculate the ratio            
            % calculate edge used ratio
            r_m_n1 =0;  % variable to hold number of mm-wave used
            r_e_n1 = 0; % variable to hold number of ethernet link used
            c_m_n1 = 0; % count total number of mm-wave link in the network topology
            c_e_n1 = 0; % count total number of ethernet link in the network topology       
            all_n1 = 0; % count total number of link
            
            mmwave_used_link_capacity =0.0; % sum of the used mm-wave link capacity
            ethernet_used_link_capacity = 0.0 ;% sum of the used ethernet link capacity
            all_used_link_capacity=0.0; % sum of the used link (mm-wave + ethernet) capacity
            
            % Calculate the mmwave edge              
            mmWavesEdges = zeros(V_expanded,V_expanded);% edges stored here
            expanded_connectivity_graph = A_con{1,1};
            for id = 1: size(ID,2)
                % find node with mmWave edge
                if ~isempty(strfind(ID{id},'_2'))
                    % find edges attached with this node and add it to mmWaves
                    mmWavesEdges(id,:) = expanded_connectivity_graph(id,:);   % provide total rows of id
                    mmWavesEdges(:,id) = expanded_connectivity_graph(:,id);    % provide total columns of id                
                end
            end
            
            normalized_mmwave_link_capacity= 0.0;
            normalized_all_link_capacity=0.0;
            normalized_eth_link_capacity = 0.0;
            
            for ii=1:V_expanded
                for jj=1:V_expanded
                    %edge is mmWave
                    if mmWavesEdges(ii,jj)==1 
                        c_m_n1 = c_m_n1 + expanded_connectivity_graph(ii,jj); % number of mm wave link  , here 'AA' matrix contains original node connections with multiple links
                        r_m_n1 = r_m_n1 + (used_edge_matrix(ii,jj)>0); % number of mm wave used 
                        all_n1 = all_n1 + (used_edge_matrix(ii,jj)>0); %number of all edges used
                        
                        % calculate the used capacity 
                        if(used_edge_matrix(ii,jj)>0)
                            normalized_mmwave_link_capacity = used_edge_matrix(ii,jj)/ max_link_capacity(ii,jj);
                            normalized_all_link_capacity = used_edge_matrix(ii,jj)/ max_link_capacity(ii,jj);
                        end
                        mmwave_used_link_capacity = mmwave_used_link_capacity + normalized_mmwave_link_capacity; % sum of used mm-wave link capacity
                        all_used_link_capacity = all_used_link_capacity + normalized_mmwave_link_capacity; % sum of used mm-wave link capacity

                    else
                        c_e_n1 = c_e_n1 + expanded_connectivity_graph(ii,jj); % number of eth link
                        r_e_n1 = r_e_n1 + (used_edge_matrix(ii,jj)>0); % number of eth link used
                        all_n1 = all_n1 + (used_edge_matrix(ii,jj)>0);
                        
                        % calculate the used capacity 
                        if(used_edge_matrix(ii,jj)>0)
                            normalized_eth_link_capacity = used_edge_matrix(ii,jj)/ max_link_capacity(ii,jj);
                            normalized_all_link_capacity = used_edge_matrix(ii,jj)/ max_link_capacity(ii,jj);
                        end
                        ethernet_used_link_capacity = ethernet_used_link_capacity + normalized_eth_link_capacity; % sum of used ethernet link capacity
                        all_used_link_capacity = all_used_link_capacity + normalized_all_link_capacity; % sum of used mm-wave link capacity
                    end
                end
            end

            %% calculate the average VNF used capacity after placing all the services for allocating n flows
             
%             sum_NFV_used_capacity=0.0;
%                normalized_NFV_used_capacity = 0.0;
%                  for n = 1 : size(M,1)
%                     for v = 1 : V_expanded
%                         for a = 1 :  size(M,3)
%                             % service n is being used at core a of node v
%                             if (served_flows_by_XPU(n,v,a) > 1)
%                                 normalized_NFV_used_capacity = (U(n,v,a) - ServiceCapacity(n,v,a))/U(n,v,a);
%                                 sum_NFV_used_capacity = sum_NFV_used_capacity+normalized_NFV_used_capacity; 
%                             end
%                         end
%                     end
%                  end
                 
            % UPDATING the result and save it 
            index = i/5;
            
%              % calculate avg delay 
%                avg_flow_allocation_delay(1,index) = avg_flow_allocation_delay(1,index)+(avg_delay/f);
%                
%                %get the cost function 
%                avg_cost(1,index) = avg_cost(1,index)+ (TotalCost);
                
               
               
            num_used_edge_after_all_simulation(1,index) = num_used_edge_after_all_simulation(1,index)+sum(sum(used_edge_matrix>0));
            percentage_of_used_mmwave_edge(1,index) = percentage_of_used_mmwave_edge(1,index)+r_m_n1;
            percentage_of_used_ethernet_edge(1,index) = percentage_of_used_ethernet_edge(1,index)+r_e_n1;
            
           % statistics of average link capacity utilization 
           if(sum(sum(used_edge_matrix>0))>0)
               linkcapacity = all_used_link_capacity/(sum(sum(used_edge_matrix>0)));
           end 
           if(r_m_n1>0)
               mmwave_capacity = mmwave_used_link_capacity/r_m_n1;
           end
           if(r_e_n1>0)
               eth_capacity = ethernet_used_link_capacity / r_e_n1;
           end
           link_capacity_utilization(1,index) = link_capacity_utilization(1,index)+linkcapacity;
           mmwave_link_capacity_utilization(1,index) = mmwave_link_capacity_utilization(1,index)+mmwave_capacity;
           ethernet_link_capacity_utilization(1,index) = ethernet_link_capacity_utilization(1,index)+eth_capacity;
                              
%             %also get the confidence interval

%           link_utilization_matrix(index,i)  = linkcapacity;
%           mmwave_link_utilization_matrix(index, i) = mmwave_capacity;
%           eth_link_utilization_matrix(index, i)= eth_capacity;


            
            total_mmwave_link=c_m_n1;
            total_eth_link=c_e_n1;

          %  Unsuccessful_flow_statistics
            unsuccessfull_flow(1,index) = unsuccessfull_flow(1,index)+number_of_unsuccessful_flow_allocation;
            
            % service placement/reuse statistics
            service_placed(1,index) = service_placed(1,index)+number_of_total_placed_services;
            service_reused(1,index) = service_reused(1,index)+number_of_total_reused_services;
            service_demand(1,index) = service_demand(1,index)+total_requested_services;
 
            % path length statistics
            path_length_avg(1,index) = path_length_avg(1,index)+path_length;
            
            % calculate avg delay 
            avg_flow_allocation_delay(1,index) = avg_flow_allocation_delay(1,index)+(avg_delay/i);
            
            %get the cost function 
            avg_cost(1,index) = avg_cost(1,index)+ (parametric_cost);
            
            %average NFV capacity utilization
            %   nfv_capacity_utilization(1,index) = sum_NFV_used_capacity/sum(sum(sum(served_flows_by_XPU)));
                
           
            % make directory to save results
            %nameDir = strcat('N',num2str(V_),'F',num2str(i));
            %mkdir(nameDir);

            %name = strcat('Output_N',num2str(V_),'_F',num2str(i),'_I',num2str(run),'.mat');
            %save(strcat(nameDir,'/',name));
            end  % if mod(i,5) == 0  .... store result after every 5 flow run
end % end of 1st for loop (traversing all flows)/ end of a single simulation
%% End of a single simulation for a random graph.store the results in temporary variables (to store finally in the global variables after finishing all the runs)
% update the global output variables 
total_service_demand = (number_of_total_placed_services + number_of_total_reused_services);

end % end of running all the simulation
disp 'time';
 toc
%% End of all the simulation for a random graph. output/average value over all the simulation run
% get the cost statistics
avg_cost_mat = avg_cost/number_of_simulation;
%percentage_of_cost = (avg_cost_mat/max_parametric_cost);
Cost_Matrix(str,:) =  avg_cost_mat;

%calculate average delay
Delay_Statistics(str,:) = avg_flow_allocation_delay/number_of_simulation;

% number of unsuccessful flow
avg_unsuccessful_flow = unsuccessfull_flow/number_of_simulation;
Unsuccessful_flow_statistics(str,:)= avg_unsuccessful_flow;  % average

% average path length
Pathlength_per_simulation= path_length_avg/number_of_simulation;
path_length_per_flow =(Pathlength_per_simulation/no_flows);  %avg_used_edge
Average_Path_Length_Global(str,:)= path_length_per_flow;  %


% percentage of service placement 
avg_percentage_of_service_placement = service_placed/number_of_simulation;
per_service_placement =(avg_percentage_of_service_placement/total_service_demand)*100;  %avg_used_edge
Service_Placement_Statistics(str,:)= per_service_placement;  %

% percentage of service reuse
avg_percentage_of_service_reuse = service_reused/number_of_simulation;
per_service_reuse =(avg_percentage_of_service_reuse/total_service_demand)*100;  %avg_used_edge
Service_Reuse_Statistics(str,:)= per_service_reuse;  %

%percentage_of_used_edge
avg_used_edge_per_simulation = num_used_edge_after_all_simulation/number_of_simulation
avg_used_edge_by_totaledge =(avg_used_edge_per_simulation/total_number_of_edge)*100;  %avg_used_edge
Global_Used_Edge_Matrix(str,:)= avg_used_edge_by_totaledge;%per_avg_used_edge

% percentage of used mm-wave link
avg_used_mmwave_edge = percentage_of_used_mmwave_edge/number_of_simulation;
per_avg_used_mmwave_edge =(avg_used_mmwave_edge/total_mmwave_link)*100;  %avg_used_edge;
Global_Used_MMWAVE_Edge_Matrix(str,:)= per_avg_used_mmwave_edge;  % average percentage of mm-wave

% percentage of used ethernet link
avg_used_eth_edge = percentage_of_used_ethernet_edge/number_of_simulation;
per_avg_used_eth_edge =(avg_used_eth_edge/total_eth_link)*100;  %avg_used_edge;
Global_Used_Ethernet_Edge_Matrix(str,:)= per_avg_used_eth_edge;  % average percentage of mm-wave

%service chain length
avg_service_demand = service_demand/number_of_simulation;
%per_service_placement =(avg_percentage_of_service_placement/total_service_demand)*100;  %avg_used_edge
Service_Demand_Statistics(str,:)= avg_service_demand;%per_service_placement;  %    

% percentage of service placement 
avg_percentage_of_service_placement = service_placed/number_of_simulation;
%per_service_placement =(avg_percentage_of_service_placement/total_service_demand)*100;  %avg_used_edge
Service_Placement_Statistics(str,:)= avg_percentage_of_service_placement;%per_service_placement;  %

% percentage of service reuse
avg_percentage_of_service_reuse = service_reused/number_of_simulation;
%per_service_reuse =(avg_percentage_of_service_reuse/total_service_demand)*100;  %avg_used_edge
Service_Reuse_Statistics(str,:)=avg_percentage_of_service_reuse;% per_service_reuse;  %

% average NFV capacity utliziation 
Final_Avg_NFV_Capacity_Utilization_Matrix(str,:) = nfv_capacity_utilization/number_of_simulation;


% average link utilization (all : mm-wave + ethernet)
avg_used_capacity_usage_per_simulation = link_capacity_utilization/number_of_simulation;
Final_Avg_Link_Capacity_Utilization_Matrix(str,:)= avg_used_capacity_usage_per_simulation;%per_avg_used_edge
   
% average link utilization (mm-wave)
avg_mmwave_used_link_capacity_per_simulation = mmwave_link_capacity_utilization/number_of_simulation;
Final_Avg_Link_Capacity_Utilization_MMWAVE_Matrix(str,:)= avg_mmwave_used_link_capacity_per_simulation;%per_avg_used_edge
    
%average link utilization (ethernet)
avg_eth_used_link_capacity_per_simulation = ethernet_link_capacity_utilization/number_of_simulation;
Final_Avg_Link_Capacity_Utilization_Ethernet_Matrix(str,:)= avg_eth_used_link_capacity_per_simulation;%per_avg_used_edge      
end  % end of the for loop running the number of graph
%% END of simulation for all the random graphs
outfilename = strcat('FinalOutput_Standard_Node',num2str(V_),'.mat');
save(outfilename, 'Cost_Matrix','Delay_Statistics','Service_Demand_Statistics','Final_Avg_NFV_Capacity_Utilization_Matrix', 'Final_Avg_Link_Capacity_Utilization_Matrix', 'Final_Avg_Link_Capacity_Utilization_MMWAVE_Matrix', 'Final_Avg_Link_Capacity_Utilization_Ethernet_Matrix','Average_Path_Length_Global', 'Service_Placement_Statistics', 'Service_Reuse_Statistics', 'Global_Used_Edge_Matrix', 'Global_Used_MMWAVE_Edge_Matrix','Global_Used_Ethernet_Edge_Matrix', 'Unsuccessful_flow_statistics');
end  % end of the function