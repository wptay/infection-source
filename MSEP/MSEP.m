function [estimated_sources, estimated_sources_count] = MSEP(G)
% G: ajacency matrix of the infection graph, i.e., the minimal tree that contains all infected nodes

%% assumptions
% assume that the maximum number of infection sources kmax = 3

%% Are there 3 sources?
% use Infection Partitioning algorithm to partition the infection graph into 3 regions
d = graphallshortestpaths(G);   % compute pariwise distance of all nodes in G
diameter = max(max(d));
sources_count_threshold = diameter/5;   % sources_count_threshold depends on the network propery of the underlying network, shall be ajusted for different kinds of networks 
n = size(G,1); % number of nodes in G
[centers, regions] = IP(G,d,3);    % use IP to partition the infection graph into 3 regions, and find the rumor center of each region

estimated_sources_count = 3;
while true
    % combine region 1 and region 2
    combined_region = union(regions{1}, regions{2});
    [Ga,Gb,w] = find(G); % Ga and Gb represents linking relationship in G
    Ga_other_region = ~ismember(Ga,combined_region);
    Ga(Ga_other_region) = [];
    Gb(Ga_other_region) = [];
    Gb_other_region = ~ismember(Gb,combined_region);
    Ga(Gb_other_region) = [];
    Gb(Gb_other_region) = [];
    G_combined_region = sparse(Ga,Gb,ones(length(Ga),1));
    % check whether the combined region is connected
    d_combined_region = graphallshortestpaths(G_combined_region);   % compute pariwise distance of all nodes in G_combined_region
    other_region = 1:1:n;
    other_region = setdiff(other_region, combined_region);
    d_combined_region(other_region,:)=0;
    d_combined_region(:,other_region)=0;
    if max(max(d_combined_region)) ~= inf
        TSE_estimates = TSE(G_combined_region, combined_region(1));
        if d(TSE_estimates(1), TSE_estimates(2)) < sources_count_threshold
            estimated_sources_count = 2;
            break;
        end
    end
        
    % combine region 2 and region 3
    combined_region = union(regions{2}, regions{3});
    [Ga,Gb,w] = find(G); % Ga and Gb represents linking relationship in G
    Ga_other_region = ~ismember(Ga,combined_region);
    Ga(Ga_other_region) = [];
    Gb(Ga_other_region) = [];
    Gb_other_region = ~ismember(Gb,combined_region);
    Ga(Gb_other_region) = [];
    Gb(Gb_other_region) = [];
    G_combined_region = sparse(Ga,Gb,ones(length(Ga),1));
    % check whether the combined region is connected
    d_combined_region = graphallshortestpaths(G_combined_region);   % compute pariwise distance of all nodes in G
    other_region = 1:1:n;
    other_region = setdiff(other_region, combined_region);
    d_combined_region(other_region,:)=0;
    d_combined_region(:,other_region)=0;
    if max(max(d_combined_region)) ~= inf
        TSE_estimates = TSE(G_combined_region, combined_region(1));
        if d(TSE_estimates(1), TSE_estimates(2)) < sources_count_threshold
            estimated_sources_count = 2;
            break;
        end
    end

    % combine region 1 and region 3
    combined_region = union(regions{1}, regions{3});
    [Ga,Gb,w] = find(G); % Ga and Gb represents linking relationship in G
    Ga_other_region = ~ismember(Ga,combined_region);
    Ga(Ga_other_region) = [];
    Gb(Ga_other_region) = [];
    Gb_other_region = ~ismember(Gb,combined_region);
    Ga(Gb_other_region) = [];
    Gb(Gb_other_region) = [];
    G_combined_region = sparse(Ga,Gb,ones(length(Ga),1));
    % check whether the combined region is connected
    d_combined_region = graphallshortestpaths(G_combined_region);   % compute pariwise distance of all nodes in G
    other_region = 1:1:n;
    other_region = setdiff(other_region, combined_region);
    d_combined_region(other_region,:)=0;
    d_combined_region(:,other_region)=0;
    if max(max(d_combined_region)) ~= inf
        TSE_estimates = TSE(G_combined_region, combined_region(1));
        if d(TSE_estimates(1), TSE_estimates(2)) < sources_count_threshold
            estimated_sources_count = 2;
            break;
        end
    end
    
    estimated_sources = centers;
    break;
end

%% Are there 2 sources?
if estimated_sources_count == 2
    [centers, regions] = IP(G,d,2);    % use IP to partition the infection graph into 2 regions, and find the rumor center of each region
    region = 1:1:n;
    TSE_estimates = TSE(G,region);
    if d(TSE_estimates(1), TSE_estimates(2)) < sources_count_threshold
        estimated_sources_count = 1;
        estimated_sources = distance_center(region, d);
    else
        estimated_sources = centers;
    end
end
