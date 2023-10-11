function distance_center = distance_center(region, d)
%   Input   d: distance matrix. d(i,j) denotes the distance between node i and node j
%           region: nodes in the considered region
%
%   Output  distance_center: set of nodes with minimum distance centrality 
%
%   Description:    distance centrality of node v = \sum_{i \in region} d(v,i)
%
n = length(d);

other_region = 1:1:n;
other_region = setdiff(other_region, region);
d(other_region,:)=0;

distance_centrality = sum(d);
distance_centers = find(distance_centrality == min(distance_centrality));
distance_center = distance_centers(ceil(rand*length(distance_centers)));


        




