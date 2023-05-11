%%
%============================================================================================================
%  Function : Check_Path_Goodness(shortestpathlist,XPU_Core_Matrix)
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
% (a) shortestPath - shortes path from a source to a destination
% (b) XPU_Core_Matrix - 
% (c) graph_link_capacity_matrix
% (d) graph_link_latency_matrix
% (e) flow_bandwidth_demand
% (f) flow_latency_limit

%============
% OUTPUT : 
%============
% 1. Integer variable 'flag_path_goodness' that can hold following two values: 
%       value 1 = Path is good. It fulfills all the constrains, good for further checking
%       value 0 = Path is not good. It does not fulfill all the constrains, no need for further checking
%
% 2.  Array 'XPU_List' of dimension (1 X Path_Length), 
%       This array holds only the 'node id' in the shortest path that are the XPUs
%============================================================================================================
%%
function [flag_path_goodness, XPU_List] = Check_Path_Goodness(shortestPath,XPU_Core_Matrix,graph_link_capacity_matrix, graph_link_latency_matrix, flow_bandwidth_demand, flow_latency_limit)
    
    P= num2str(shortestPath); % converting from number vector to string
    Path = str2num(P); % converting the string 
    %disp 'considering path';
    
    %% defining variables
    flag_path_goodness=0;
    No_of_XPU =0;      
    Flag_link_capacity_latency =1; 
	
	XPU_List = []; %zeros(1,length(Path)); %%
	
    
    start_index=1; % not considering the source node
    next_index = (start_index + 1);
    end_index = length(Path); % not considering the destination node 
    
    %% traverser the shortest path to check if this path fulfills all the constrain
    while (start_index <= end_index)  
         % check if the node is an XPU
         if(XPU_Core_Matrix(shortestPath(1,start_index))>0)
            No_of_XPU=No_of_XPU+1;
			XPU_List(1, No_of_XPU) = shortestPath(1,start_index);  %%
         end
        
         if(start_index ~= end_index)
             
             current_node = Path(start_index);
             next_node = Path(next_index);
             % check link capacity        
             if(graph_link_capacity_matrix(current_node,next_node)<flow_bandwidth_demand)
                 Flag_link_capacity_latency=0;
                 break;
             end
             
             % check link latency
             if(graph_link_latency_matrix(current_node,next_node)>flow_latency_limit)
                 Flag_link_capacity_latency=0;
                 break;
             end
         end %end of if(start_index ~= end_index)
         
         %increase index
         start_index = start_index+1; 
         next_index = (start_index + 1);
    end  % end of the while loop  
    
    % path goodness checking - final decision for this path
    if((No_of_XPU>0) && (Flag_link_capacity_latency ==1))       
        flag_path_goodness =1;
        %{
         else
        if(No_of_XPU == 0)
            disp 'NO XPU';
        end
        if(Flag_link_capacity_latency ==0)
            disp 'not available link bandwidth / low latency';
        end
        %}
    end    
end % end of the function
    
         