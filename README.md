# ResourceManage5G
%% HOW TO RUN THE CODE:

For multiple run on multi-graphs
--------------------------------------

Configuration of the script
=============================
The script that should be run is 'runnerFunction_RMA_Heuristic_MultiGraph_latest.m'. Before running the script we need to put the input files under folder 'Input'. 

Input :   
========
The input files should be kept in folder named 'Input'. The pre-generated input files for node 10, 15, 20, and 25 are stored under folder 'Multigraph-Input'. 
For instance, if we need to run the simulation for 30 nodes for 10 different graphs, we should search for the input files for 30 nodes in folder 'Multigraph-Input' and copy those to folder 'Input'.  


Output
=======
The output files are stored in the same location where the script file is located.  

Run the Script
====================
In the Matlab console the script 'runnerFunction_RMA_Heuristic_MultiGraph_latest.m' should run with the following parameters:

runnerFunction_RMA_Heuristic_MultiGraph_latest( number of nodes, number of simulation to run)

for instance, to run 50 simulations with 30 nodes (for 10 different topologies) and from 1 to 100 flows to allocate, from Matlab prompt just paste the following line and wait:
runnerFunction_RMA_Heuristic_MultiGraph_latest( 30, 50)
