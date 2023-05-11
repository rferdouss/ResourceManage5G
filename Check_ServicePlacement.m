%%
%============================================================================================================
%  Function : Check_ServicePlacement(pathcost, path_XPU_List, Service_chain, Max_core_per_XPU, NFV_Bandwidth_Matrix,flow_service_demand, service_placement_cost_matrix, service_matrix, service_allocation_xpu, service_allocation_xpu_core);
%  Date :  August 2016
%  -----------------------------------------------------------------------------
%% This function checks the possibility of satisfying the flow-graph (service chain demain) along a shortest path 
% This function takes a shortest path as input and considering the current network status check the possibility of satisfying the flow-graph (service chain demain) along a shortest path with counting the reuse of services or placing new services in XPU along the path. 
% Lets assume, the shortest path is sp, flow f demand for the service chain (4->2->3->1), two cases may occur:
%(1) Reuse of service is possible :
%(2) All new services should be placed: 


%============
% INPUT parameters: 
%============
% (a) pathcost - cost of a shortest path
% (b) path_XPU_List - list of XPUs in this shortest path
% (c) Service_chain - service chain demand of this flow to be served
% (d) Max_core_per_XPU - Maximum number of core avaialable pe XPU
% (e) NFV_Bandwidth_Matrix - 3d Matrix listing the available bandwidth after placing some NFV/services
% (f) flow_service_demand - Demand of the flow for services
% (g) service_placement_cost_matrix - Demand of the flow for services
% (h) service_matrix - Demand of the flow for services
% (i) service_allocation_xpu - Demand of the flow for services
% (j) service_allocation_xpu_core - Demand of the flow for services

%============
% OUTPUT : 
%============  , , , 
% 1. Integer variable 'flag_flowallocation' that can hold following two values: 
%       value 1 = Path is good for fullfilling the service chain demand. 
%       value 0 = Path is not good. It does not fulfill all the service chain demand, can be discard
%
% 2.  'totalpathcost' containing the updated path cost after fullfilling the service chain/flow graph
%
% 3.  array 'service_allocation_xpu_core' containing XPUs core id order that can fulfill the flow graph.
%     for example, if the service chain demand is 3->1->4 and the XPU core ids can be 1->2->3
% 4.  array 'service_allocation_xpu' containing XPUs id order that can fulfill the flow graph.
%     for example, if the service chain demand is 3->1->4 and the XPU ids can be a->b->b
%============================================================================================================
%%
%% FUNCTION
function [flag_flowallocation, totalpathcost, service_allocation_xpu, service_allocation_xpu_core] = Check_ServicePlacement(pathcost, path_XPU_List, Service_chain, Max_core_per_XPU, NFV_Bandwidth_Matrix,flow_service_demand, service_placement_cost_matrix, service_matrix, service_allocation_xpu, service_allocation_xpu_core);
%% Defining variables
number_of_placed_service =0; % count the number of placed services
totalpathcost=0;  % calculate the total cost of flow allocation along this shortest path (link dynamic cost + link fixed cost + service placement cost)
 %% Case 2 : 
  %CODE CHECKING -  if the remaining service is possible to palce in this path
    if (any(service_allocation_xpu)>0) % if atleast one service can be reused
    %    disp 'reuse is possible';
     %   disp 'current service allocation is as below';
      %  service_allocation_xpu
       for i=1: length(service_allocation_xpu)
       %    disp 'starting loop --status'
        %   i
         %  service_allocation_xpu
          % disp 'now considering';
           %service_allocation_xpu(1,i)
           %disp 'number of placed services'
          % number_of_placed_service
           if(service_allocation_xpu(1,i)==0)
           %    disp 'it is 0 so here';
               
              % define variables
              next_xpu_node_id = 1000;
              next_xpu_core = 1000;
              prev_xpu_node_id=0;
              prev_xpu_core=0;
               
              subset = service_allocation_xpu(i:length(service_allocation_xpu)); 
              subset_xpu_core = service_allocation_xpu_core(i:length(service_allocation_xpu_core)); 
              %indices = find(subset, 1, 'first') % find the indices insided the array containing the first non-zero value   
              %indices_core =  find(subset_xpu_core, 1, 'first') % find the indices insided the array containing the first non-zero value 
              
              if(isempty(find(subset, 1, 'first'))==0)
                  next_xpu_node = subset(1, find(subset, 1, 'first'));
                  next_xpu_node_id  = find(path_XPU_List == next_xpu_node);
                  next_xpu_core = subset_xpu_core(1, find(subset_xpu_core, 1, 'first'));
              end
              if(i==1)
                  prev_xpu_node=0;
                  prev_xpu_node_id=0;
                  prev_xpu_core=0;
              else
                  prev_xpu_node = service_allocation_xpu(1,(i-1));
                  prev_xpu_node_id = find(path_XPU_List == prev_xpu_node);
                  prev_xpu_core = service_allocation_xpu_core(1,(i-1));
              end
              % find a suitable xpu node to place this service
              for j=1 : length(path_XPU_List)
                  node_placed=0;
                  % check if the XPU order supports placing this services
                %  if((path_XPU_List(1,j)<=next_xpu_node) && (next_xpu_node>0))
                  if((j<=next_xpu_node_id) && (next_xpu_node_id>0))
                      
                       if((j>=prev_xpu_node_id))
                      %if((path_XPU_List(1,j)>=prev_xpu_node))
                          for  k=1:Max_core_per_XPU % checking the cores of XPU
                              % if core id is the one where the previous or next service is already placed
                              contains_xpu_core = any(service_allocation_xpu_core(:) == k);
                             
                              if(contains_xpu_core==0)
                                  if (service_matrix(Service_chain(i, 1),path_XPU_List(1,j),k)==1)
                                    %  service_matrix(Service_chain(i, 1),path_XPU_List(1,j),k)
                                      % check if enough bandwidth is available 
                                      if(NFV_Bandwidth_Matrix(Service_chain(i, 1),path_XPU_List(1,j), k)>flow_service_demand)
                                        %  NFV_Bandwidth_Matrix(Service_chain(i, 1),path_XPU_List(1,j), k)
                                         % disp 'service placement is possible in this XPU node core';
                                          totalpathcost = totalpathcost + service_placement_cost_matrix(Service_chain(i, 1),path_XPU_List(1,j), k);
                                          service_allocation_xpu(1, i)= path_XPU_List(1,j);                                         
                                        
                                          service_allocation_xpu_core(1, i)= k;
                                        %  disp 'before adding a new services';
                                         % number_of_placed_service
                                          number_of_placed_service=number_of_placed_service+1;
                                        %  disp 'a new service is placed';
                                         % number_of_placed_service
                                          node_placed=1;
                                          break;
                                          %flag_serviceplaced=1;                     
                                      end
                                  end
                              end
                          end
                      end
                  end
                  if(node_placed==1)
                      break;
                  end
              end % end of  for j=1 :  length(path_XPU_List)
           else
               number_of_placed_service=number_of_placed_service+1;
             %  disp 'else if service is already placed';
             %  number_of_placed_service
           end % end if(service_allocation_xpu(1,i)==0) 
           
          % disp 'number of placed service at the end';
          % number_of_placed_service
       end  % end of the for loop   
    %   disp 'after running the loop the service allocation becomes';
     %  service_allocation_xpu
    end
    
%% Case 1 :  currently no service is placed in this path, check if the required flow-demand (service chain) is possible to place in this path 
 if (all(~service_allocation_xpu)) 
   %  disp 'need all service placement';
     start_xpu_index=1;  % start index of the XPU list
     start_xpu_core = 1; % start index of the XPU core
     flag_serviceplaced=0;

     for i = 1: length(Service_chain)
         while (start_xpu_index <= length(path_XPU_List)) % traverse all the XPUs
             while (start_xpu_core <= Max_core_per_XPU)
                 % check if the service i is possible to place in this core
                 if(service_matrix(Service_chain(i, 1),path_XPU_List(1,start_xpu_index), start_xpu_core)==1) 
                     % check if enough bandwidth is available 
                     if(NFV_Bandwidth_Matrix(Service_chain(i, 1),path_XPU_List(1,start_xpu_index), start_xpu_core)>flow_service_demand)
                       %  disp 'service placement is possible in this XPU node core';
                         totalpathcost = totalpathcost + service_placement_cost_matrix(Service_chain(i, 1),path_XPU_List(1,start_xpu_index), start_xpu_core);
                         service_allocation_xpu(1, i)= path_XPU_List(1,start_xpu_index);
                    
                         service_allocation_xpu_core(1, i)= start_xpu_core;
                         number_of_placed_service=number_of_placed_service+1;
                         flag_serviceplaced=1;                     
                     end
                 end
                 
                 if (start_xpu_core == Max_core_per_XPU) % increase the xpu
                     start_xpu_index = start_xpu_index+1;
                     start_xpu_core =1;
                     break;
                 else
                     start_xpu_core = start_xpu_core+1;
                 end
                 
                 if(flag_serviceplaced==1)
                     break;
                 end
             end  % end of the inner while loop
             
             % if a service/NFG is fulfilled/reused then break the loop right here and go to see if any other next service can be reused
             if(flag_serviceplaced==1)
                 flag_serviceplaced=0;
                 break;
             end
         end  % end of the outer while loop 
     end % end of the for loop 
 end


 %% Decision on flow allocation and service placement  - calculate the total path cost
     totalpathcost = totalpathcost+pathcost; 
   %  disp 'number of service placed/reuse';
   %  number_of_placed_service
     
     if(number_of_placed_service == length(Service_chain))
         flag_flowallocation =1;
    %     if(any(service_allocation_xpu(1,:)==0))
     %        disp '';
      %   end
  
     else
         flag_flowallocation=0;
     end
end
 