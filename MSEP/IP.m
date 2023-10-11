function [centers, regions] = IP(G,d, roots_count)

%% inputs
% G: ajacency matrix of the infection graph, i.e., the minimal tree that contains all infected nodes
% d: pairwise distance of all nodes in G
% roots_count: number of regions to partition into

%% initial guess of the centers
% first find one diameter
diameter = max(max(d));
[start_leaf, end_leaf] = find(d==diameter);
diameter = spath(G, start_leaf(1), end_leaf(1));
% then find the intial guess of centers 
switch roots_count
    case 2
        centers=[diameter(round((length(diameter)+1)/3)) diameter(round((length(diameter)+1)/3*2))]; % initial guess: 1/3 and 2/3 of the diameter
    case 3
        centers=[diameter(round((length(diameter)+1)/4)) diameter(round((length(diameter)+1)/4*2)) diameter(round((length(diameter)+1)/4*3))]; % initial guess: 1/4, 2/4 and 3/4 of the diameter
end

%% iteration start
end_iteration = false;
error_tolerence = 1;
iteration_count = 0;
while end_iteration == false
    new_centers = zeros(1,roots_count);
    distance_between_new_and_current_centers = zeros(1,roots_count);
    distance_to_centers = d(centers,:);
    minimum_distance_to_centers = min(distance_to_centers);
    regions = cell(1,roots_count);
    for i = 1:roots_count
        regions{i} = find(distance_to_centers(i,:) == minimum_distance_to_centers);  % partition G into roots_count regions based on current centers
        new_centers(i) = distance_center(regions{i}, d);  % find the distance center of each region (for trees, distance center and rumor center has been proved to be the same)
        distance_between_new_and_current_centers(i) = d(new_centers(i), centers(i));
    end
    centers = new_centers;
    iteration_count=iteration_count+1;
    if max(distance_between_new_and_current_centers) <= error_tolerence || iteration_count >= 20
        end_iteration = true;
    end
end
end

        
       
