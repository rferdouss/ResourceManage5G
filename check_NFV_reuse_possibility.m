%%
%=================================================================================================================================================
%  Function : check_NFV_reuse_possibility(path_XPU_List, Service_chain, Max_core_per_XPU, NFV_in_XPU, NFV_Bandwidth_Matrix, flow_service_demand)
%  Date :  August 2016
%  -----------------------------------------------------------------------------
%% This function checks the goodness of a shortest path based on availability of XPUs, enough link bandwidth and low link latency. 
% To check the goodness of the shortest path the following three constrains should be fulfilled:
% 1. Constrain 1 - there should exist at least one XPU in the way of the shortest path, i.e., at least one node of the shortest path should be an XPU
% 2. Constrain 2 - availability of enough link bandwidth throughout the shortest path than the flow demand
% 3. Constrain 3 - availability of low link latency throughout the shortest path than the flow latency
%
%============
% INPUT : 
%============
% (c) path_XPU_List - list of XPUs in this shortest path
% (d) Service_chain - service chain demand of this flow to be served
% (e) Max_core_per_XPU - Maximum number of core avaialable pe XPU
% (f) NFV_in_XPU - 3d Matrix listing the placed NFV/services
% (g) NFV_Bandwidth_Matrix - 3d Matrix listing the available bandwidth after placing some NFV/services
% (h) flow_service_demand - Demand of the flow for services
%============
% OUTPUT : 
%============
% 1. Integer variable 'flag_path_goodness' that can hold following two values: 
%       value 1 = Path is good. It fulfills all the constrains, good for further checking
%       value 0 = Path is not bad. It does not fulfill all the constrains, no need for further checking
%
% 2.  Array 'XPU_List' of dimension (1 X Path_Length), 
%       This array holds only the 'node id' in the shortest path that are the XPUs

% This function takes a shortest path as input, checks the possibility of satisfying the service (NFV) flowgraph. 
% Output :
% a) On success - 
%		i) flag =1, indicating that this path is a potential candidate for allocating this flow
%		ii) cost>0 - total cost of using this path and flow allocation 
% a) On faulure - 
%		i) flag = 0, indicating that this path is not a potential candidate for allocating this flow
%		ii) cost = 0, 

%  Output  : 
% Two flags : 1. flag_NFV_Reuse_possibility = flag to indicate if the
% service placement is possible or not
% 2. longest number of service reuse =
% 3. longest sub-set of the services can be reused


% UPDATE :  This function should return should chech 
% 1. Number of highest possible services that can be reused
% 2. If the other services can be possible to placed in this path
% 3. Final output is the flag (1/0), and the total cost of the path

%=========================================================================================================================================================
%%  [flag_flow_allocation, totalcost, number_of_service_reuse, service_allocation_xpu, service_allocation_xpu_core
function [flag_NFV_Reuse_possibility, totalcost, longest_service_chain_length, longest_service_chain, longest_service_chain_xpu_core] = check_NFV_reuse_possibility(pathcost, path_XPU_List, Service_chain, Max_core_per_XPU, NFV_in_XPU, NFV_Bandwidth_Matrix, flow_service_demand, service_placement_cost_matrix, service_matrix)
    
%% defining variables
    
    
    
    %% Find the largest sub-set among the service chain that can be reused (for example, for a flow f the service demand is 4->2->3->1. Lets assume, the longest service chain that can be fulfilled is 3->1)
    
    flag_NFV_Reuse_possibility= 0;  
    totalcost=0;
      
    longest_service_chain = zeros(1, length(Service_chain));  % array to hold the longest sub-set of the service chain that can be reused
    longest_service_chain_xpu_core = zeros(1, length(Service_chain));  % array to hold the longest sub-set of the service chain that can be reused
    longest_service_chain_length=0; % number of the services in the longest sub-set of the service chain that can be reused
   
    
    for i=1 :length(Service_chain)
        current_service_chain_allocation = zeros(1, length(Service_chain));  % array to hold the longest sub-set of the service chain that can be reused
        current_service_chain_allocation_xpu_core = zeros(1, length(Service_chain));
        current_service_chain_length=0; % number of the services in the longest sub-set of the service chain that can be reused  
        current_xpu_id=1;
        previous_xpu_id=0;   
        
        for service_chain_id =i : length(Service_chain)
            for j =1 : length(path_XPU_List)
                current_xpu_id = j;%path_XPU_List(1,j);  % consider the current xpu
                % check all the cores of this XPU to find if this service is already placed there and has enough bandwidth to serve the new flow
                for k =1 : Max_core_per_XPU
                    % check if the service i is available in this core
                    if(NFV_in_XPU(Service_chain(service_chain_id, 1),path_XPU_List(1,j), k)==1)
                        % check if enough bandwidth is available                         
                        if(NFV_Bandwidth_Matrix(Service_chain(service_chain_id, 1),path_XPU_List(1,j), k)>flow_service_demand)
                            %disp 'service is already placed here and it is possible to use this service';
                            if (current_xpu_id >=previous_xpu_id)
                                current_service_chain_allocation(1, service_chain_id)= path_XPU_List(1,j);
                                current_service_chain_allocation_xpu_core(1, service_chain_id)= k;
                                current_service_chain_length = current_service_chain_length+1;
                                previous_xpu_id = current_xpu_id;
                                break;
                            end
                        end
                    end
                end % end of the most- inner for loop - for j =1 : Max_core_per_XPU
            end % end of the inner for loop
        end     % end of for service_chain_id =i : length(Service_chain) 
        
       % disp 'before sending to function check_serviceplacement';
       % current_service_chain_allocation
        [flag_flowallocation, totalpathcost, service_allocation_xpu, service_allocation_xpu_core] = Check_ServicePlacement(pathcost, path_XPU_List, Service_chain, Max_core_per_XPU, NFV_Bandwidth_Matrix,flow_service_demand, service_placement_cost_matrix, service_matrix, current_service_chain_allocation, current_service_chain_allocation_xpu_core);
       % disp 'return flow allocation value';
        %flag_flowallocation
        %UPDATE - Monday 05/09/2016 - a path will be consider the longest
        %service chain if the remaining services can be placed there,
        %otherwise discard
        if (flag_flowallocation ==1)
          
         %    if(any(service_allocation_xpu(1,:)==0))
          %   disp 'service allocation returns with 0';
         %end
            
            if(i==1)
                longest_service_chain = service_allocation_xpu;
                longest_service_chain_xpu_core = service_allocation_xpu_core;
                longest_service_chain_length = current_service_chain_length;
                totalcost = totalpathcost;
                flag_NFV_Reuse_possibility = flag_flowallocation;
                %check if the other service is possible to place here
            else
                if (current_service_chain_length>longest_service_chain_length)
                    % if multiple longest service chains of same length are found then take the one with minimal residual bandwidth
                  longest_service_chain = service_allocation_xpu;
                  longest_service_chain_xpu_core = service_allocation_xpu_core;
                  longest_service_chain_length = current_service_chain_length;  
                  totalcost = totalpathcost;
                  flag_NFV_Reuse_possibility = flag_flowallocation;
                end
            end
        end
        %UPDATE : Check for the total cost after REUSING ANY SERVICE AND PLACING NEW SERVICE - this function should return, 1.flag_possibility of reusing any service and placement, 2. total cost, number of service to be reused, final service placement xpu and core matrix
        % TOOK the minimum cost solution and at the same time the flow-graph (service chain) should be fulfilled
    end % end of the outer for loop
    
   % flag_NFV_Reuse_possibility = flag_flowallocation;
    
    %DECISION ABOUT 'flag_NFV_Reuse_possibility' IS TAKEN AFTER CHECKING IF
    %THE REMAINING services can be placed here or not
    %CODE : 
    % X = [1 0 4 8 0 0 0 8 6];
    % X(2:5)
    % k = X(4:length(X))
    % indices = find(k, 1, 'first')
   
    %   number_of_existing_service_use=0; % count the number of NFV/Services can be reused for this flow
  %  service_chain_id =1;
    
    %%
    % algorithm description : Friday 26/08/2016
    % Input : A shortest path with at least one XPU
    % Process : 
    %       1. Find the longest service chain (among the total service chain) that is possible to reuse 
    %       2. if size(longest service chain)  == size (service chain demand)
    %                       flow is served in this path
    %       3. check other possibilities with placing other remaining services and re-using all or a set of already placed services
%% check if the already placed NFV/services can be reused for satisfying the current flow service chain demand
%{
 for i =1 : length(path_XPU_List)
        while (service_chain_id <=length(Service_chain))
            for j =1 : Max_core_per_XPU
                % check if the service i is available in this core
                if(NFV_in_XPU(Service_chain(1, service_chain_id),path_XPU_List(1,i), j)==1)
                    % check if enough bandwidth is available 
                    if(NFV_Bandwidth_Matrix(Service_chain(1, service_chain_id),path_XPU_List(1,i), j)>flow_service_demand)
                        disp 'service is already placed here and it is possible to use this service';
                        number_of_existing_service_use=number_of_existing_service_use+1;
                        service_chain_id =service_chain_id+1; % go for the next service to be served
                        break;
                    end
                end
            end   % end of the inner for loop 
            if(j==Max_core_per_XPU)  % if all the cores in an XPU has been examined, then break it and go for the next XPU 
                break;
            end 
        end % end of the while loop
    end % end of the for loop
    
    if(number_of_existing_service_use == length(Service_chain))
        flag_NFV_Reuse_possibility=1;
    end
end % end of the function
   
%}
    
    %% check if the already placed NFV/services can be reused for satisfying the current flow service chain demand
    %{
     start_xpu_index=1;  % start index of the XPU list
    start_xpu_core = 1; % start index of the XPU core
   
       for i = 1: length(Service_chain)
        while (start_xpu_index <= length(path_XPU_List)) % traverse all the XPUs
            while (start_xpu_core <= Max_core_per_XPU)
                % check if the service i is available in this core
                if(NFV_in_XPU(Service_chain(1, i),path_XPU_List(1,start_xpu_index), start_xpu_core)==1)
                    % check if enough bandwidth is available 
                    if(NFV_Bandwidth_Matrix(Service_chain(1, i),path_XPU_List(1,start_xpu_index), start_xpu_core)>flow_service_demand)
                        disp 'service is already placed here and it is possible to use this service';
                        number_of_existing_service_use=number_of_existing_service_use+1;
                        break;
                    end
                end
                
                if (start_xpu_core == Max_core_per_XPU) % increase the xpu
                    start_xpu_index = start_xpu_index+1;
                    start_xpu_core =1;
                    break;
                else
                    start_xpu_core = start_xpu_core+1;
                end
            end  % end of the inner while loop
            
            % if a service/NFG is fulfilled/reused then break the loop right here and go to see if any other next service can be reused
            if(number_of_existing_service_use==i)
                break;
            end           
        end  % end of the outer while loop   
    end  % end of the for loop
    %}
end     
  
        


    
         