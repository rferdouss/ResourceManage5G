function [ A_, X, Y ] = functionGenerateGraph_leo_v2(N,d,T,range_mm, area_x, area_y, cp, l, c_c, c_d)

C_ethernet = 2*10^8;
C_mmWave = 3*10^8;

% link cost/capcity etc
% capacity
cp_eth = cp(1,1);
cp_mm = cp(1,2);
cp_o = cp(1,3);

% latency
l_eth = l(1,1);
l_mm = l(1,2);
l_o = l(1,3);

% cost fixed
c_c_eth = c_c(1,1);
c_c_mm = c_c(1,2);
c_c_o = c_c(1,3);

% cost dynamic
c_d_eth = c_d(1,1);
c_d_mm = c_d(1,2);
c_d_o = c_d(1,3);


% no of technologies used
no_tech = size(T,2);

% adjuncy matrix
A = zeros(N,N,no_tech);
% fixed code
A_c = zeros(N,N,no_tech);
% dynamic cost
A_d = zeros(N,N,no_tech);
% capacity
A_cp = zeros(N,N,no_tech);
% latency
A_l = zeros(N,N,no_tech);

%% FORMULATION

% total number of links
S = (d*N)/2;

% check if it is possible to make a graph
if (S < N-1) || (S > N*(N-1)/2) || T(1,1) < 0 || T(1,1) > 1 || T(1,2) < 0 || T(1,2) > 1 || T(1,3) < 0 || T(1,3) > 1
    A = -1;
    return;
end

coordinates=zeros(N,2);

% put points
for i=1:N
    x = rand()*area_x;
    y = rand()*area_y;
    coordinates(i,:)= [x,y];
end

% put mmwave edges between the nodes which are close
count_mm = 0;
count_eth = 0;

for i = 1:N
    for j = 1:N
        
        if i == j
            continue;
        end
        
        % find distance between node i and j
        dis = sqrt( (coordinates(i,1)-coordinates(j,1))^2 + (coordinates(i,2)-coordinates(j,2))^2 );
        
        % can have an mmWave connection
        if dis <= range_mm
            A(i,j,1)=1; % put ethernet
            A(i,j,2)=1; % put mmWave
            
            A_c(i,j,1)= c_c_eth;
            A_c(i,j,2)= c_c_mm ;
            
            A_d(i,j,1)=c_d_eth;
            ps = successProbV3(dis,'3gpp', 'no');
            A_d(i,j,2)= 1/ps;
            
            A_cp(i,j,1)=cp_eth;
            A_cp(i,j,2)=cp_mm;
                        % tranmission delay    + propagation delay
            A_l(i,j,1)= (1500*8)/(cp_eth*10^6) + dis/C_ethernet;
                        % tranmission delay    + propagation delay
            A_l(i,j,2)= (1500*8)/(cp_mm*10^6) + dis/C_mmWave;
            %A_l(i,j,2)= dis/C_mmWave;
            
            count_eth=count_eth+0.5;
            count_mm=count_mm+0.5;
        end
        
    end
end

G = graph(A(:,:,1));
bin = conncomp(G);

while max(bin) ~= 1
    for ii = 1:max(bin)
        
        ind_s = find(bin==ii);
        
        %%%
        % connect nodes which are closest in these two clusters
        min = [-1 -1];
        min_dis = inf;
        
        
        for jj=1:max(bin)
            if ii==jj
                continue;
            end
            ind_t = find(bin==jj);
            
            for a = 1: size(ind_s,2)
                for b = 1: size(ind_t,2)
                    i = ind_s(a);
                    j = ind_t(b);
                    
                    i_x = coordinates(i,1);
                    i_y = coordinates(i,2);
                    
                    j_x = coordinates(j,1);
                    j_y = coordinates(j,2);
                    
                    dis = sqrt( (i_x-j_x)^2 + (i_y-j_y)^2 );
                    
                    if dis <= min_dis
                        min_dis = dis;
                        min = [i j];
                    end
                end
            end
        end
        i = min(1,1);
        j = min(1,2);
        dis = min_dis;
        %%%
        
        % out an eth edge between then
        A(i,j,1)=1; % put ethernet
        A(j,i,1)=1; % put ethernet
        
        A_c(i,j,1)= c_c_eth;
        A_c(j,i,1)= c_c_eth;
        
        A_d(i,j,1)=c_d_eth;
        A_d(j,i,1)=c_d_eth;
        
        A_cp(i,j,1)=cp_eth;
        A_cp(j,i,1)=cp_eth;
        
        A_l(i,j,1)= (1500*8)/(cp_eth*10^6) + dis/C_ethernet;
        A_l(j,i,1)= (1500*8)/(cp_eth*10^6) + dis/C_ethernet;
        
        count_eth=count_eth+1;
    end
    G = graph(A(:,:,1));
    bin = conncomp(G);
end

count_eth = sum(sum(A(:,:,1)))/2;

%% randomly connect nodes now in the graph with ethernet edges
no_edges = N;
for a = 1 : N
    
    while true
        i = randi(N);
        j = randi(N);
        if i ~= j
            break;
        end
    end
    
    dis = sqrt( (coordinates(i,1)-coordinates(j,1))^2 + (coordinates(i,2)-coordinates(j,2))^2 );
    
    % out an eth edge between then
    A(i,j,1)=1; % put ethernet
    A(j,i,1)=1; % put ethernet
    
    A_c(i,j,1)= c_c_eth;
    A_c(j,i,1)= c_c_eth;
    
    A_d(i,j,1)=c_d_eth;
    A_d(j,i,1)=c_d_eth;
    
    A_cp(i,j,1)=cp_eth;
    A_cp(j,i,1)=cp_eth;
    
    A_l(i,j,1)=(1500*8)/(cp_eth*10^6) + dis/C_ethernet;
    A_l(j,i,1)=(1500*8)/(cp_eth*10^6) + dis/C_ethernet;
    
end

%G = graph(A(:,:,1));
%bin = conncomp(G);

X = coordinates(:,1);
Y = coordinates(:,2);
%[G h] = drawGraph( A , X, Y, area_x, area_y);

A_ = {A,A_c,A_d,A_cp,A_l };


end

