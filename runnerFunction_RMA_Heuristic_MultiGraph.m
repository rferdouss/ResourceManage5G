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
function [  ] = runnerFunction_RMA_Heuristic_MultiGraph( V_, number_of_simulation)

number_of_graph =10;
no_flows = 100; 

%% Global matrix to hold the final result (over m run on f flows on n number of graphs)
Global_Used_Edge_Matrix= zeros(number_of_graph,20);  % PERCENGAGE  of link (mm-wave + ethernet) utilization
Global_Used_MMWAVE_Edge_Matrix= zeros(number_of_graph,20);  % percentage of mm-wave link utilization
Global_Used_Ethernet_Edge_Matrix= zeros(number_of_graph,20);  % percentage of ethernet link utilization

Unsuccessful_flow_statistics = zeros(number_of_graph,20);  % number of unsuccessfull flow
Average_Path_Length_Global = zeros(number_of_graph,20);  % average path length
% service demand
Service_Placement_Statistics = zeros(number_of_graph,20);  % Number of total placed services
Service_Reuse_Statistics = zeros(number_of_graph,20);  % number of total reused services

for str =1 : number_of_graph
    disp '----------------------------------------------------------------------------';
    fprintf('Random Graph ID : %d\n',str);
    %disp  'Next RANDOM GRAPH';
    disp '-----------------------------------------------------------------------------';

%% declaration of input variable
%disp('Number of nodes:')
%disp(V_)


%% LOAD INPUT FILE (with pre-generated flow)
V = V_;             % Number of nodes e.g., [10 15 20 25 30]
    % Number of flows
filename ='';
filename = strcat('Input/Input_N_',num2str(V),'_F_',num2str(no_flows),'_i_',num2str(str),'.mat'); % input file with the generated flow, network toptology, capacity, latency, service allocation information
%filename
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
F_video = FlowGenerator(V_,no_flows_video,M_size, d, d_sd, l, l_sd);

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
F_voice = FlowGenerator(V_,no_flows_voice,M_size, d, d_sd, l, l_sd);

F = [F_video, F_voice];
% random permutation of flows- Shuffle flows
F = {F{randperm(length(F))}};

%% Make M where services can be places! Given the total number of XPU and XFE nodes
%% Metric O, which is number of cores available on each node
V_expanded = size(A_con{1,1},1);

%disp('Expanded nodes :')
%disp(V_expanded) 

%disp 'Total number of links/edges : ';
%sum(sum(A_con{1, 1}))

ratio_XPU = 0.3;
% no of XPU nodes
no_XPU = floor(V_*ratio_XPU);
% no of XFE nodes
no_XFE = V_ - no_XPU;
% no of cores available on XPU nodes
no_core_XPU = 4;

M = zeros(M_size,V_expanded, no_core_XPU);

O = zeros(1,V_expanded);

for i=1:no_XPU
    while true
        % select a node
        r = randi(V_);
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



%% Semi-Globar Variables (that remain active till the end of all simulation run for a single topology/graph), e.g., 50 runs for 100 flows with n nodes  
k =30;                                  % NUMBER OF shortest path to investigate from source node to destination
num_used_edge_after_all_simulation= zeros(1,20);
percentage_of_used_mmwave_edge = zeros(1,20);
percentage_of_used_ethernet_edge = zeros(1,20);


service_reused = zeros(1,20);
service_placed = zeros(1,20);
total_service_demand=0;

total_mmwave_link=0;
total_eth_link=0;

% 
unsuccessfull_flow = zeros(1,20);

path_length_avg = zeros(1,20);

%% RUNNING FOR MULTIPLE INSTANCES OF ALLOCATING ALL THE FLOWS AMONG THE NODES : from srt to end instances
for run = 1:number_of_simulation
 
 %fprintf('Simulation run ID : %d \n', run);
%% Input Graph from loaded input file
y = A_con{1,1}; % graph connectivity

total_number_of_edge = sum(sum(y));
%total_number_of_edge
% fixed cost of using link
y_c =A_con{1,2};
% dynamic cost of using link
y_d =A_con{1,3};
% link capacity
y_cp =A_con{1,4};
graph_capacity = A_con{1,4};
% link latency
y_l =A_con{1,5};
% no of nodes
N = size(y,1); % Number of nodes

%disp 'Number of nodes after graph expansion';
%N
%% INPUT Flow parameters

ServiceMatrix =  M;
ServiceCapacity = U;
ServiceCostMatrix = h;


%% Defining other local variable (active throughout all the flow allocation) for a single simulation run
NFV_Placed_in_XPU = zeros(M_size,V, no_core_XPU);  % keep track of all the placed NFVs/Services
used_edge_matrix=zeros(V,V); % matrix to store the used edge information 

% statistics about NFV placement and reuse
number_of_total_placed_services = 0;
number_of_total_reused_services = 0;
path_length=0;%path length
parametric_cost=0;%parametric cost
number_of_unsuccessful_flow_allocation=0;

%% Initial preparation before running the loop of flow allocation 
%first prepare the cost matrix (in order to run the K-Shortest path algorithm)
y_temp= y_c+y_d;  % summing up the dynamic and fixed cost 
y_temp( y_temp==0 )=Inf;  % for the requirement of the k-shortest path calculation code all the cost that is zero are replaced by 'INF', indicating the absence of coxt/edge
%y_temp
% Randomization of the flow by shuffle all flows
F = {F{randperm(length(F))}};
   
%% RUNNING MAIN LOOP - ALLOCATION OF ALL FLOW AND SERVICE PLACEMENT -- FOR THIS SIMULATION RUN
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
   
    %% Defining the variables  (local). works for a single flow
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
  %  fprintf('The path request is from source node %d to destination node %d, with K = %d \n',S,D, k);
 % disp 'shortes path number';
%	 length(shortestPaths) 
     
 %    if( length(shortestPaths) <30)
  %       disp 'flow id';
   %      i
    %     disp 'simulation id';
     %    run
      %   disp 'graph id';
       %  str
         
       %   [shortestPaths, totalCosts] = kShortestPath(y_temp, S, D, k);
  %   end
     
	%%    
    %Flow allocation - running over the k shortest path to find the optimum one with placing required/desired services
    if isempty(shortestPaths)
      %  fprintf('No shortest path available between these nodes\n\n');
    else
        for j = 1: length(shortestPaths)  % 2nd for loop (running over all the k shortest path to find the optimum one for allocation this flow)
        %    fprintf('\n'); 
         %  fprintf('Path # %d:\n',j);
          % disp(shortestPaths{j})
           % disp (S_chain);
          %  fprintf('Cost of path %d is %5.2f\n\n',j,totalCosts(j));
           
            %Defining local variables
            service_allocation_among_xpu = zeros (1,length(S_chain)); % keep the id of the xpu from the array 'XPU_List' that contains the list of XPUs
            service_allocation_among_xpu_core = zeros (1,length(S_chain));
            number_of_service_reuse=0;
            flag_service_placement_possibility=0; % this flag indicates if the services can be placed along this shortest path or not
          %  path_cost = 0;
            
            %% FILTER LEVEL 1 : CHECK GOODNESS OF THE SHORTEST PATH (Constrain 1: exists at least one XPU, and Constrain 2. availability of enough link bandwidth and low latency)
            
            [Flag_Path_Goodness, XPU_List] = Check_Path_Goodness(shortestPaths{j},O,graph_capacity, y_l, d_f, l_f);
           
          %  disp 'xpu exist result';
           % Flag_Path_Goodness
            %% FILTER LEVEL 2: CHECK THE POSSIBILITY OF FLOW ALLOCATION IN THE PATH - if this shortest path is good (contains at least one XPU, enough link bandwidth, and low latency) then go for next level to check if this shortest path has the potential to fulfil the flow-graph
           
            if(Flag_Path_Goodness==1)
               % disp 'Filter 1 : Path is good. It contains at least one XPU. Have sufficient bandwidth and low latency- going to check more for the possibility of satisfying the flow-graph ';   %[flag_Optimalpath, path_cost, service_allocation_xpu, service_allocation_xpu_core] = Check_Flow_Allocation_Possibility(shortestPaths{j}, totalCosts(j), XPU_List, S_chain, no_core_XPU, NFV_Placed_in_XPU, U, d_f, ServiceCostMatrix, ServiceMatrix); %check the possibility of service placement and flow allocation in this shortest path         
                if (all(~NFV_Placed_in_XPU))   % if this is the first flow/ no service is placed yet/ no existing service  
                   % disp 'No NFV is placed yet / NFV_in_XPU matrix is empty';    
                    [flag_flowallocation, totalcost,service_allocation_xpu, service_allocation_xpu_core]= Check_ServicePlacement(totalCosts(j), XPU_List, S_chain, no_core_XPU, ServiceCapacity,d_f, ServiceCostMatrix, ServiceMatrix, service_allocation_among_xpu, service_allocation_among_xpu_core);
                    service_allocation_among_xpu = service_allocation_xpu;
                    service_allocation_among_xpu_core = service_allocation_xpu_core;
                    
                  %  disp 'result from othe function';
                   % service_allocation_xpu;
                    %service_allocation_xpu_core
                    %disp 'result local';
                    %service_allocation_among_xpu
                    %service_allocation_among_xpu_core
                    flag_service_placement_possibility = flag_flowallocation;
                   % path_cost = totalcost;
                   % fprintf('%d \n', flag_flowallocation);
                else % next flow/ NFV_in_XPU matrix is not empty / go for checking if existing services can be reused   
                    % Check 1 : Check for best option - No new service placement - Use of peviously placed services   
                    [flag_service_placement_possibility, totalcost, number_of_service_reuse, service_allocation_xpu, service_allocation_xpu_core] = check_NFV_reuse_possibility(totalCosts(j),XPU_List, S_chain, no_core_XPU, NFV_Placed_in_XPU, ServiceCapacity, d_f,ServiceCostMatrix, ServiceMatrix); %check the possibility of service placement and flow allocation in this shortest path         
                end  %end of if (all(~NFV_in_XPU)) 
                
                %% Now check the result of the service placement possibility -  if this path is a candidate for optimal path
                 if(flag_service_placement_possibility==1) % service placement is possible in this path
                 %   disp 'Filter 1 : this shortest path is a candidate to be an optimal path for this flow';
                    % for the first shortest path, assigns it to optimal cost and path for checking with the next shortest paths
                    %disp 'result from othe function';
                    %service_allocation_xpu;
                    %service_allocation_xpu_core
                    %disp 'result local';
                    %service_allocation_among_xpu
                    %service_allocation_among_xpu_core
                    
                    
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
                        
                        
                       %optimal_shortestPaths
                        %optimal_service_allocation_among_xpu
                        %optimal_service_allocation_among_xpu_cores
                    end
                    
                    flag_Optimalpath=1;
                 else %  if(flag_service_placement_possibility==0)
                    % COMMENT: disp 'Filter 1 : Service placement error - All Service placement is not possible in this path';
                 end 
      
            else
            %COMMENT    disp 'Filter 1: Service placement is not possible in this path- No XPU node/Low link bandwidth/High latency';
				%flag_Optimalpath=0;
            end 
          %  disp 'path id';
          %  j
            
           % disp 'length shortest path';
            %length(shortestPaths);
        end % end of 2nd for loop (traversing the set of shortest paths found for a flow from source to destination node)				
       
        
        %disp 'optimal path result';
        %flag_Optimalpath
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
                  path_length = path_length+length(optimal_shortestPaths);
                  parametric_cost =parametric_cost+optimal_Costs; 
                  %used_edge_matrix
                  
        else
          %  fprintf('FLOW NO : %d,  START NODE %d, END NODE %d \n', i, S, D);
          %COMMENT  disp 'Filter 2 : This flow allocation is not possible - no optimal path found among the list of K shortest paths';
          %COMMENT  disp 'This flow allocation is not possible'; 
          disp '---------';
          fprintf('flow id = %d, simulation id =%d, graph id = %d\n', i, run, str);
          fprintf('source %d, destination %d, total shortest path %d \n',S,D, length(shortestPaths));
            number_of_unsuccessful_flow_allocation=number_of_unsuccessful_flow_allocation+1;
           % y
            
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
            
            for ii=1:V_expanded
                for jj=1:V_expanded
                    %edge is mmWave
                    if mmWavesEdges(ii,jj)==1 
                        c_m_n1 = c_m_n1 + expanded_connectivity_graph(ii,jj); % number of mm wave link  , here 'AA' matrix contains original node connections with multiple links
                        r_m_n1 = r_m_n1 + (used_edge_matrix(ii,jj)>0); % number of mm wave used 
                        all_n1 = all_n1 + (used_edge_matrix(ii,jj)>0); %number of all edges used
                    else
                        c_e_n1 = c_e_n1 + expanded_connectivity_graph(ii,jj); % number of eth link
                        r_e_n1 = r_e_n1 + (used_edge_matrix(ii,jj)>0); % number of eth link used
                        all_n1 = all_n1 + (used_edge_matrix(ii,jj)>0);
                    end
                end
            end

            % UPDATING the result and save it 
            index = i/5;
            
            num_used_edge_after_all_simulation(1,index) = num_used_edge_after_all_simulation(1,index)+sum(sum(used_edge_matrix>0));
            percentage_of_used_mmwave_edge(1,index) = percentage_of_used_mmwave_edge(1,index)+r_m_n1;
            percentage_of_used_ethernet_edge(1,index) = percentage_of_used_ethernet_edge(1,index)+r_e_n1;
            
            total_mmwave_link=c_m_n1;
            total_eth_link=c_e_n1;

          %  Unsuccessful_flow_statistics
            unsuccessfull_flow(1,index) = unsuccessfull_flow(1,index)+number_of_unsuccessful_flow_allocation;
            
            % service placement/reuse statistics
            service_placed(1,index) = service_placed(1,index)+number_of_total_placed_services;
            service_reused(1,index) = service_reused(1,index)+number_of_total_reused_services;
 
            % path length statistics
            path_length_avg(1,index) = path_length_avg(1,index)+path_length;
           
             % make directory to save results
          %  nameDir = strcat('N',num2str(V_),'F',num2str(i));
           % mkdir(nameDir);

            %name = strcat('Output_N',num2str(V_),'_F',num2str(i),'_I',num2str(run),'.mat');
            %save(strcat(nameDir,'/',name));
            end  % if mod(i,5) == 0  .... store result after every 5 flow run
end % end of 1st for loop (traversing all flows)/ end of a single simulation
%% End of a single simulation for a random graph.store the results in temporary variables (to store finally in the global variables after finishing all the runs)
% update the global output variables 
total_service_demand = (number_of_total_placed_services + number_of_total_reused_services);

end % end of running all the simulation
%% End of all the simulation for a random graph. output/average value over all the simulation run

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


end  % end of the for loop running the number of graph
%% END of simulation for all the random graphs

outfilename = strcat('FinalOutput_Node',num2str(V_),'.mat');
save(outfilename, 'Average_Path_Length_Global', 'Service_Placement_Statistics', 'Service_Reuse_Statistics', 'Global_Used_Edge_Matrix', 'Global_Used_MMWAVE_Edge_Matrix','Global_Used_Ethernet_Edge_Matrix', 'Unsuccessful_flow_statistics');
end  % end of the function



