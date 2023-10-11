function distance_center = function_distance_center(explicit, d)
%   Input   d: distance matrix. d(i,j) denotes the distance between node i and node j
%           explicit: array of explicit nodes
%
%   Output  distance_center: set of nodes with minimum distance centrality 
%
%   Description:    distance centrality of node v = \sum_{i \in V_e} d(v,i)
%
n = length(d);

non_explicit = 1:1:n;
non_explicit = setdiff(non_explicit, explicit);
d(non_explicit,:)=0;

distance_centrality = sum(d);
distance_center = find(distance_centrality == min(distance_centrality));


        




