function function_find_shortest_path(neighbors, d, pre, s, t, L)
global shortest_paths 
if d(s,t) == 1
    shortest_paths(end+1,:) = [pre s t];
end

if d(s,t) > 1 
    all_neighbors = neighbors{s};
    for counter = 1:length(all_neighbors)
        neighbor = all_neighbors(counter);
        if d(neighbor, t) < L
            function_find_shortest_path(neighbors, d, [pre s], neighbor, t, L-1);
        end
    end
end