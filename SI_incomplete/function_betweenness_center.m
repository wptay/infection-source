function betweenness_center = function_betweenness_center(explicit, d, neighbors)
global shortest_paths 
n = length(d);
betweenness_centrality = zeros(1,n);
for i=1:length(explicit)-1
    s = explicit(i);
    for j=i+1:length(explicit)
        t = explicit(j);
        shortest_paths = [];
        pre = [];
        function_find_shortest_path(neighbors, d, pre, s, t, d(s,t));
        k_set = unique(shortest_paths);
        for counter=length(k_set)
            k = k_set(counter);
            if k ~= s && k~= t
                betweenness_centrality(k) = betweenness_centrality(k)+length(find(shortest_paths == k))/size(shortest_paths,1);
            end
        end
    end
end
betweenness_center = find(betweenness_centrality == max(betweenness_centrality));

        




