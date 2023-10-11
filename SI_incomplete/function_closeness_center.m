function closeness_center = function_closeness_center(explicit, d)
%   Input   d: distance matrix. d(i,j) denotes the distance between node i and node j
%           explicit: array of explicit nodes
%
%   Output  closeness_center: set of nodes with maximum closeness centrality 
%
%   Description:    closeness centrality of node v = \sum_{i \in V_e, i \ne v} \frac{1}{d(v,i)}
%
n = length(d);
d(logical(eye(size(d)))) = Inf;     %set the diagonal elements to be infinity, later when we take inverse, it becomes 0

non_explicit = 1:1:n;
non_explicit = setdiff(non_explicit, explicit);
d(non_explicit,:)=Inf;

d = 1./d;
closeness_centrality = sum(d);
closeness_center = find(closeness_centrality == max(closeness_centrality));


        




