%% This script generates the input files that can used for the simulation

clear all
clc
close all
 
% {
%% Generate Graph
plotG = 'yes';
% no of nodes
V_org = 30;
% number of flows
no_flows = 100;


% node degree
d = 3; % NOT USED
% percentage of links
T = [1 0.5 0]; % NOT USED

% link cost/capcity etc
% capacity [ethernet mmWave optical]
cp = [10000 2000 0];% Mb/s

% latency [ethernet mmWave optical]
l = [1 1 0]; % dependent on distance= d/c

% cost fixed [ethernet mmWave optical]
c_c = [100 1 0];

% cost dynamic [ethernet mmWave optical]
c_d = [1 1 0]; % NOT USED for mmWave! Calculated inside the function

% range distance for mmWave link
range_mm = 200;

% area we are intertested in
area_x = 2000;
area_y = 2000;

%% Generate the Initial Graph
[A, xx, yy] = functionGenerateGraph_leo_v2(V_org,d,T,range_mm, area_x, area_y, cp, l, c_c, c_d);

%% Generate ids of the node. You can add them manually as well
ID = {};
for i = 1:V_org
    name = strcat(num2str(i));
    ID = [ID name];
end

%% Convert Graph: so that each two nodes in the graph share only one link atmost
ID_or = ID;
[A_con, ID] = convertGraph(A, ID);
V = size(A_con{1,1},1);

%% plot graphs to see results
% {

if plotG == 'yes'
    
    % Original Graph generated
    AA = A{1};
    drawGraph( AA , xx, yy, area_x, area_y);
    
    % color converted graph
    n_color = zeros(V,3);
    n_color(:,3) = 1;
    for i = 1 : V
        % mark virtual nodes
        if i > V_org
            n_color(i,:) =[0 1 0];
        end
    end
    
    % original graph with multi links
    A_adjuncy = A{1,1};
    AA = zeros(size(A_adjuncy,1),size(A_adjuncy,1));
    for i=1:size(A_adjuncy,1)
        for j=1:size(A_adjuncy,1)
            AA(i,j) = sum(A_adjuncy(i,j,:));
        end
    end
    
    G = graph(AA, ID_or);
    weights = G.Edges.Weight;
    e_color_vector = zeros(size(weights,1),3);
    for e = 1:size(weights,1)
        if weights(e) > 1
            e_color_vector(e,:) = [1 0 0];
        else
            e_color_vector(e,:) = [0 0 1];
        end
    end
    
    % graph with multi-links
    figure
    p = plot(G,'EdgeLabel',G.Edges.Weight);
    p.EdgeColor = e_color_vector;
    
    % graph without multi-links
    G_con = graph(A_con{1,1}, ID);
    figure
    p = plot(G_con);
    p.NodeColor = n_color;
    
end

%% Generate Network Flows

% ratio of flows that are video flows
videoFlows = 0.70;

% flows for video and voice
no_flows_video = round(no_flows*0.7);
no_flows_voice = no_flows-no_flows_video;

% number of services in the network
M_size = 4;

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
F_video = FlowGenerator(V_org,no_flows_video,M_size, d, d_sd, l, l_sd);

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
F_voice = FlowGenerator(V_org,no_flows_voice,M_size, d, d_sd, l, l_sd);

F = [F_video, F_voice];
% random permutation of flows- Shuffle flows
F = {F{randperm(length(F))}};

%% Make M where services can be places! Given the total number of XPU and XFE nodes
%% Metric O, which is number of cores available on each node
ratio_XPU = 0.3;
% no of XPU nodes
no_XPU = floor(V_org*ratio_XPU);
% no of XFE nodes
no_XFE = V_org - no_XPU;
% no of cores available on XPU nodes
no_core_XPU = 4;

M = zeros(M_size,V, no_core_XPU);

O = zeros(1,V);

for i=1:no_XPU
    while true
        % select a node
        r = randi(V_org);
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

U = zeros(M_size,V, no_core_XPU);
for i = 1:M_size
    U(i,:,:) = no_flow_eachVNF;
end

%% Create h(n,i), fixed cost of placing VNF of type n on node i
cost = 1000;
h = zeros(M_size, V, no_core_XPU) + cost;


%% Save the input file
filename = strcat('Input_N_',num2str(V_org),'_F_',num2str(no_flows),'.mat');
save(filename);
V
V_org







