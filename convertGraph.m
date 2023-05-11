function [ AA, ID ] = convertGraph(A_, ID)
%convertGraph Converts 3D graph A_ into 2D graphs.
% Input: A = {A,A_c,A_d,A_cp,A_l}, all are 3D array
% Output:AA= {A_con,A_c_con,A_d_con,A_cp_con,A_l_con}

%% Sample Input
% A_ = {A,A_c,A_d,A_cp,A_l};
A = A_{1,1};
A_c = A_{1,2};
A_d = A_{1,3};
A_cp = A_{1,4};
A_l = A_{1,5};

V = size(A,1);

A_con = zeros(V,V);
A_c_con = zeros(V,V);
A_d_con = zeros(V,V);
A_cp_con = zeros(V,V);
A_l_con = zeros(V,V);

for i = 1 : V
    for j = 1 : i-1
        if sum(A(i,j,:)) > 1
            
            for t = 1 : 3
                if A(i,j,t) ~= 1
                    continue;
                end
                
                % add the new node name
                node_name = strcat(ID(i),ID(j),'_',num2str(t));
                ID = [ID, node_name];
                
                % convert the main array to 2D array
                A_con_new = zeros(size(A_con,1)+1,size(A_con,1)+1);
                A_con_new(1:size(A_con,1),1:size(A_con,1)) = A_con;
                A_con_new(i,j) = 0;
                A_con_new(size(A_con,1)+1,i) = 1;
                A_con_new(size(A_con,1)+1,j) = 1;
                A_con = A_con_new;
                
                % convert the fixed cost array to 2D array
                A_c_con_new = zeros(size(A_c_con,1)+1,size(A_c_con,1)+1);
                A_c_con_new(1:size(A_c_con,1),1:size(A_c_con,1)) = A_c_con;
                A_c_con_new(i,j) = 0;
                % fixed cost is divided between the links
                A_c_con_new(size(A_c_con,1)+1,i) = A_c(i,j,t)/2;
                A_c_con_new(size(A_c_con,1)+1,j) = A_c(i,j,t)/2;
                A_c_con = A_c_con_new;
                
                % convert the dynamic cost array to 2D array
                A_d_con_new = zeros(size(A_d_con,1)+1,size(A_d_con,1)+1);
                A_d_con_new(1:size(A_d_con,1),1:size(A_d_con,1)) = A_d_con;
                A_d_con_new(i,j) = 0;
                % dynamic cost is divided between the links
                A_d_con_new(size(A_d_con,1)+1,i) = A_d(i,j,t)/2;
                A_d_con_new(size(A_d_con,1)+1,j) = A_d(i,j,t)/2;
                A_d_con = A_d_con_new;
                
                % convert the capacity array to 2D array
                A_cp_con_new = zeros(size(A_cp_con,1)+1,size(A_cp_con,1)+1);
                A_cp_con_new(1:size(A_cp_con,1),1:size(A_cp_con,1)) = A_cp_con;
                A_cp_con_new(i,j) = 0;
                % capacity is same on both links
                A_cp_con_new(size(A_cp_con,1)+1,i) = A_cp(i,j,t);
                A_cp_con_new(size(A_cp_con,1)+1,j) = A_cp(i,j,t);
                A_cp_con = A_cp_con_new;
                
                % convert the latency array to 2D array
                A_l_con_new = zeros(size(A_l_con,1)+1,size(A_l_con,1)+1);
                A_l_con_new(1:size(A_l_con,1),1:size(A_l_con,1)) = A_l_con;
                A_l_con_new(i,j) = 0;
                % latency is divided between the links
                A_l_con_new(size(A_l_con,1)+1,i) = A_l(i,j,t)/2;
                A_l_con_new(size(A_l_con,1)+1,j) = A_l(i,j,t)/2;
                A_l_con = A_l_con_new;
                
            end
        else
            % its either zero or one, so you add it to the final graph
            A_con(i,j) = sum(A(i,j,:));
            % update the link costs and everything
            A_c_con(i,j) = A_c(i,j);
            A_d_con(i,j) = A_d(i,j);
            A_cp_con(i,j) = A_cp(i,j);
            A_l_con(i,j) = A_l(i,j);
        end
    end
end

for i = 1 : size(A_con,1)
    for j = 1 : i-1
        A_con(j,i) = A_con(i,j); 
        A_c_con(j,i) = A_c_con(i,j); 
        A_d_con(j,i) = A_d_con(i,j); 
        A_cp_con(j,i) = A_cp_con(i,j);
        A_l_con(j,i) = A_l_con(i,j);
    end
end

%G = graph(A_con, ID);
%figure
%h = plot(G);

AA = {A_con, A_c_con, A_d_con, A_cp_con, A_l_con };

end

