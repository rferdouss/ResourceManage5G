%% This script generates the input files that can used for the simulation

clear all
clc
close all
 
%% Generating multiple graphs 
number_of_random_graph =10;
for str=1:number_of_random_graph
% {
%% Generate Graph
plotG = 'yes';

% no of nodes
V_org = 40;
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

if (plotG == 'yes')
    
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
else
    disp '';
    
end


%% Save the input file
filename = strcat('Input/Input_N_',num2str(V_org),'_F_',num2str(no_flows),'_i_',num2str(str),'.mat');
save(filename);
disp 'Graph id : ';
str
disp 'input node';
V_org
disp 'node after expansion';
V
disp 'Total number of links/edges : ';
sum(sum(A_con{1, 1}))

end % end of the for loop






