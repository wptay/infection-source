clear;
rand('state',sum(100*clock));
%network_type = 'facebook';
%network_type = 'power_grid';
%network_type = 'small_world';
network_type = 'as';
number_of_runs=1000;
M = 200;                        %maximum number of infected nodes
MAX = 10000;                      %upper bound of infected nodes
switch network_type
    case 'small_world'          %total number of underlying nodes
        n = 1000;
    case 'power_grid'
        n = 4941;               %total number of underlying nodes in power grid network
    case 'facebook'
        n = 4039;
    case 'as'
        n = 3933;
end

explicit_percentage =0.5;      %percentage of explicit nodes wrt infected nodes, value 0 means q \in [max(0,2-1/p),1]
save_file_name=strcat('SI_Incomplete_RG_BC_CC_DC_',network_type,'_M',num2str(M),'_n',num2str(n),'_explicit_percentage',num2str(explicit_percentage*100),'_run',num2str(number_of_runs));

%all_run_G = cell(number_of_runs,1);
all_run_source = zeros(number_of_runs,1);
all_run_explicit = cell(number_of_runs,1);
all_run_infected = cell(number_of_runs,1);
all_run_diameter = zeros(number_of_runs,1);
all_run_candidate_size = zeros(number_of_runs,1);
%FOr Reverse Greedy
all_run_fval_Reverse_Greedy = cell(number_of_runs,1);
all_run_estimate_Reverse_Greedy = zeros(number_of_runs,1);
all_run_error_Reverse_Greedy = zeros(number_of_runs,1);
all_run_time_Reverse_Greedy = zeros(number_of_runs,1);
%For betweeness center
all_run_estimate_betweeness = zeros(number_of_runs,1);
all_run_error_betweeness = zeros(number_of_runs,1);
all_run_time_betweeness = zeros(number_of_runs,1);
%For closeness center
all_run_estimate_closeness = zeros(number_of_runs,1);
all_run_error_closeness = zeros(number_of_runs,1);
all_run_time_closeness = zeros(number_of_runs,1);
%For distance center
all_run_estimate_distance = zeros(number_of_runs,1);
all_run_error_distance = zeros(number_of_runs,1);
all_run_time_distance = zeros(number_of_runs,1);

if (strcmpi(network_type,'power_grid'))
    load(['datasets/power_grid/' 'power_grid_network.mat'],'G', 'd_power_grid', 'neighbors');
end
if (strcmpi(network_type,'facebook'))
    load(['datasets/facebook/' 'facebook_network.mat'],'G', 'd_facebook', 'neighbors');
end
if (strcmpi(network_type,'as'))
    load(['datasets/as/' 'as_network.mat'],'G', 'd_as', 'neighbors');
end
for counter=1:number_of_runs
    if mod(counter,5)==0
        disp(counter);
    end
    %% Generate the underlying network
    switch network_type
        case 'small_world'          %total number of underlying nodes
            G = smallw(n,1,0.3);        %small_world_network
            source=1;                   %specify the source node
            [G_a, G_b, w]=find(G);
            w = ones(1,length(G_a));
            G = sparse(G_a, G_b, w);
            neighbors = cell(n,1);
            for i=1:n
                neighbors{i} = G_b(G_a == i)';
            end
        case 'power_grid'
            source = ceil(rand*n);
        case 'facebook'
            source = ceil(rand*n);
        case 'as'
            source = ceil(rand*n);
    end
    all_run_source(counter) = source;
    
    %% Spread the infection on the underlying network
    explicit=[];
    while isempty(explicit)
        [explicit,infected]=function_infection_spread(G,M,MAX,source,explicit_percentage);
    end
    if (strcmpi(network_type,'small_world'))
        all_run_G{counter,1} = G;
    end
    all_run_explicit{counter,1} = explicit;
    all_run_infected{counter,1} = infected;
    
    %% Compute distance between any pair of nodes
    switch network_type
        case 'small_world'          
            d=[];
            for i=1:n
                [d(i,:), path, pred] = graphshortestpath(G, i);       
            end
        case 'power_grid'
            d = d_power_grid;
        case 'facebook'
            d = d_facebook;
        case 'as'
            d = d_as;
    end
    
    max_distances = max(d);
    diameter = max(max_distances);
    all_run_diameter(counter)=diameter;
    
    %% betweenness centralilty based estimator
    estimates_betweeness = function_betweenness_center(explicit, d, neighbors);
    index=max(1,ceil(length(estimates_betweeness)*rand));
    estimate_betweeness=estimates_betweeness(index);
    all_run_estimate_betweeness(counter)=estimate_betweeness;
    [all_run_error_betweeness(counter), path, pred] = graphshortestpath(G, source, estimate_betweeness);
    
    %% closeness centralilty based estimator
    estimates_closeness = function_closeness_center(explicit, d);
    index=max(1,ceil(length(estimates_closeness)*rand));
    estimate_closeness=estimates_closeness(index);
    all_run_estimate_closeness(counter)=estimate_closeness;
    [all_run_error_closeness(counter), path, pred] = graphshortestpath(G, source, estimate_closeness);
    
    %% distance centralilty based estimator
    estimates_distance = function_distance_center(explicit, d);
    index=max(1,ceil(length(estimates_distance)*rand));
    estimate_distance=estimates_distance(index);
    all_run_estimate_distance(counter)=estimate_distance;
    [all_run_error_distance(counter), path, pred] = graphshortestpath(G, source, estimate_distance);
    
    %% Reverse Greedy Algorithm 
    candidate_relax = 1;
    non_explicit = 1:1:n;
    non_explicit = setdiff(non_explicit, explicit);
    d(non_explicit,:)=0;
    max_min_distances = max(d);
    estimates=[];
    min_d=min(max_min_distances);
    for i=1:n
        if max_min_distances(i) <= min_d+candidate_relax
            estimates(end+1)=i;
        end
    end
    candidate_size = length(estimates);
    all_run_candidate_size(counter)=candidate_size;
    if candidate_size > 1
        each_run_fval_MIQCQP = zeros(n,1);
        each_run_fval_Reverse_Greedy = zeros(n,1);
        each_run_x = cell(n,1);
        for index_root=1:candidate_size
            root = estimates(index_root);
            %% Minimum connected subgraph containing root and explicit nodes
            [dist, path, pred] = graphshortestpath(G, root);
            T_sp_a = [];
            T_sp_b = [];
            for i=1:n
                j = pred(i);
                if j ~= 0 && ~isnan(j)
                    T_sp_a = [T_sp_a i j];
                    T_sp_b = [T_sp_b j i];
                end
            end
            w = ones(length(T_sp_a),1);
            T_sp = sparse(T_sp_a, T_sp_b, w);
            one_vector = ones(n,1) ;
            non_explicit_leaf = 0;
            explicit_root = [explicit root];
            G_sub = G;          %minimum connected subgraph containing root and explicit nodes
            [G_sub_a, G_sub_b, w]=find(G_sub);
            while ~isempty(non_explicit_leaf)
                one_vector = ones(size(T_sp,1),1) ;
                degree = T_sp*one_vector;
                leaf = find(degree==1);
                non_explicit_leaf = setdiff(leaf,intersect(leaf, explicit_root));
                temp = ismember(T_sp_a,non_explicit_leaf);
                T_sp_a = T_sp_a(~temp);
                T_sp_b = T_sp_b(~temp);
                temp = ismember(T_sp_b,non_explicit_leaf);
                T_sp_a = T_sp_a(~temp);
                T_sp_b = T_sp_b(~temp);
                w = ones(length(T_sp_a),1);
                T_sp = sparse(T_sp_a, T_sp_b, w);
                temp = ismember(G_sub_a,non_explicit_leaf);
                G_sub_a = G_sub_a(~temp);
                G_sub_b = G_sub_b(~temp);
                temp = ismember(G_sub_b,non_explicit_leaf);
                G_sub_a = G_sub_a(~temp);
                G_sub_b = G_sub_b(~temp);
                w = ones(length(G_sub_a),1);
                G_sub = sparse(G_sub_a, G_sub_b, w);
            end
            [head,tail,w] = find(G_sub);    %sparse representation, edge: from tail to head
            G_sub_nodes = unique(head);
           
            %% Reverse Greedy Algorithm
            tic;
            each_run_fval_Reverse_Greedy(root) = function_Reverse_Greedy(G_sub, T_sp, root);
            all_run_time_Reverse_Greedy(counter)=all_run_time_Reverse_Greedy(counter)+toc;
        end
        %For Reverse Greedy
        temp = setdiff(each_run_fval_Reverse_Greedy,0);
        min_fval = min(temp);
        estimates_Reverse_Greedy = find(each_run_fval_Reverse_Greedy==min_fval);
        index=max(1,ceil(length(estimates_Reverse_Greedy)*rand));
        estimate_Reverse_Greedy=estimates_Reverse_Greedy(index);
        all_run_fval_Reverse_Greedy{counter,1}=each_run_fval_Reverse_Greedy;
    else
        estimate_Reverse_Greedy = estimates;
    end
    all_run_estimate_Reverse_Greedy(counter)=estimate_Reverse_Greedy;
    [all_run_error_Reverse_Greedy(counter), path, pred] = graphshortestpath(G, source, estimate_Reverse_Greedy);
end

%% Plot the error distances cdf
max_error_Reverse_Greedy = max(all_run_error_Reverse_Greedy);
max_error_betweeness = max(all_run_error_betweeness);
max_error_closeness = max(all_run_error_closeness);
max_error_distance = max(all_run_error_distance);
max_error = max([max_error_Reverse_Greedy max_error_betweeness max_error_closeness max_error_distance]);
for i=1:max_error+1
    error_bar(i,1)=length(find(all_run_error_Reverse_Greedy<=i-1))/number_of_runs*100;
    error_bar(i,2)=length(find(all_run_error_betweeness<=i-1))/number_of_runs*100;
    error_bar(i,3)=length(find(all_run_error_closeness<=i-1))/number_of_runs*100;
    error_bar(i,4)=length(find(all_run_error_distance<=i-1))/number_of_runs*100;
end 
x=[0:1:max_error];
scrsz = get(0,'ScreenSize');
figure('Position',[0 30 scrsz(3) scrsz(4)-95]);
h=bar(x,error_bar,'grouped');
legend 'RG' 'BC' 'CC' 'DC';
xlabel('error distance cdf'); ylabel('percentage (%)'); 
set(h(1),'facecolor','red') % use color name
set(h(2),'facecolor',[0 1 0]) % or use RGB triple
set(h(3),'facecolor',[0 0 1]) % or use RGB triple
set(h(4),'facecolor',[0 0 0]) % or use RGB triple
saveas(gcf,[save_file_name '_cdf'],'bmp');
close(gcf);

%% Plot the error distances pdf
for i=1:max_error+1
    error_bar(i,1)=length(find(all_run_error_Reverse_Greedy==i-1))/number_of_runs*100;
    error_bar(i,2)=length(find(all_run_error_betweeness==i-1))/number_of_runs*100;
    error_bar(i,3)=length(find(all_run_error_closeness==i-1))/number_of_runs*100;
    error_bar(i,4)=length(find(all_run_error_distance==i-1))/number_of_runs*100;
end 
x=[0:1:max_error];
scrsz = get(0,'ScreenSize');
figure('Position',[0 30 scrsz(3) scrsz(4)-95]);
h=bar(x,error_bar,'grouped');
legend 'RG' 'BC' 'CC' 'DC';
xlabel('error distance pdf'); ylabel('percentage (%)'); 
set(h(1),'facecolor','red') % use color name
set(h(2),'facecolor',[0 1 0]) % or use RGB triple
set(h(3),'facecolor',[0 0 1]) % or use RGB triple
set(h(4),'facecolor',[0 0 0]) % or use RGB triple
saveas(gcf,[save_file_name '_pdf'],'bmp');
close(gcf);

%% save workspace
save([save_file_name '.mat']);