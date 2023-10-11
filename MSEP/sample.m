clear;
infection_graph = 'infection_graph\geometric_tree_sample.mat';
load(infection_graph);
%% infection graph contains the following variables
% G: ajacency matrix of the infection graph, i.e., the minimal tree that contains all infected nodes
% sources: the real infection sources 

[estimated_sources, estimated_sources_count] = MSEP(G);
