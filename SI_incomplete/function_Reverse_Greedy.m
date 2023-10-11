function fval=function_Reverse_Greedy(G, T, root)
%G is the underlying minimum connected subgraph
%T is a shortest-path tree that spans G
T=sparse(T);

% %% test graph
% clear
% G = zeros(10,10);
% G(1,2)=1;
% G(1,3)=1;
% G(1,4)=1;
% G(2,1)=1;
% G(2,3)=1;
% G(2,5)=1;
% G(2,6)=1;
% G(3,1)=1;
% G(3,2)=1;
% G(3,5)=1;
% G(3,7)=1;
% G(3,8)=1;
% G(4,1)=1;
% G(4,8)=1;
% G(5,2)=1;
% G(5,3)=1;
% G(5,6)=1;
% G(5,9)=1;
% G(5,10)=1;
% G(6,2)=1;
% G(6,5)=1;
% G(7,3)=1;
% G(8,3)=1;
% G(8,4)=1;
% G(9,5)=1;
% G(9,10)=1;
% G(10,5)=1;
% G(10,9)=1;
% T = zeros(10,10);
% T(1,2)=1;
% T(1,3)=1;
% T(1,4)=1;
% T(2,1)=1;
% T(2,5)=1;
% T(2,6)=1;
% T(3,1)=1;
% T(3,7)=1;
% T(3,8)=1;
% T(4,1)=1;
% T(5,2)=1;
% T(5,9)=1;
% T(5,10)=1;
% T(6,2)=1;
% T(7,3)=1;
% T(8,3)=1;
% T(9,5)=1;
% T(10,5)=1;
% G=sparse(G);
% T=sparse(T);
% root = 1;
% %% CLOSE Test Graph

[head,tail,w] = find(G);
G = sparse(head,tail,w);
[head_T,tail_T,w] = find(T); 
T = sparse(head_T,tail_T,w);
node_number = length(unique(head));       %number of nodes in G             
%% Compute upward-distance U 
[U, path, pa] = graphshortestpath(T, root);
%% Compute downward-distance D
n = size(G,1);
one_vector = ones(n,1);
D = zeros(1,n);
visited = 0;
L = 0;
T_temp = T;
[T_temp_a, T_temp_b, w]=find(T_temp);
while visited<node_number-1
    one_vector_temp = ones(size(T_temp,1),1) ;
    degree = T_temp*one_vector_temp;
    leaf = find(degree==1);
    leaf = setdiff(leaf,intersect(leaf,root));
    D(leaf')=L;
    temp = ismember(T_temp_a,leaf);
    T_temp_a = T_temp_a(~temp);
    T_temp_b = T_temp_b(~temp);
    temp = ismember(T_temp_b,leaf);
    T_temp_a = T_temp_a(~temp);
    T_temp_b = T_temp_b(~temp);
    w = ones(length(T_temp_a),1);
    T_temp = sparse(T_temp_a, T_temp_b, w);
    visited = visited+length(leaf);
    L = L+1;
end
D(root)=L;
%% ---------------- Get the optimal tree T_star ---------------- %%
head_star = [];
tail_star = [];
visited = [];
for level = L:-1:1
    nodes = find(U==level);     % Consider nodes at level "level" 
    %% Sort x
    if length(nodes)>1
        nodes_matrix = zeros(length(nodes),2);          %first col are the candidate nodes, 2nd col is the constraint 
        nodes_matrix(:,1)=nodes';
        nodes_matrix(:,2)=D(nodes)';
        nodes_matrix = sortrows(nodes_matrix,2);        %in order of increasing Dx
        nodes = nodes_matrix(:,1)';
    end
    %% Consider each x
    for i=1:length(nodes)
        x = nodes(i);
        visited(end+1) = x;
        %% Get candidates, and compute the D_prime value for each candidate
        Nx = tail(find(head==x));   %neighbor of x in G
        Nx = setdiff(Nx, intersect(Nx,visited));
        remove_index = [];
        D_prime = zeros(length(Nx),1);
        for j=1:length(Nx)
            y = Nx(j);
            if D(x)+U(y)+1>L
                remove_index(end+1)=j;      %index of inqualified neighbors  
            else
                if pa(x) ~= y                                 %if y is not parent of x in T
                    D_prime(j) = D(y);
                else                                            %if y is parent of x in T
                    Ny = tail_T(find(head_T==y));               %neighbor of y in T
                    Ny = setdiff(Ny,[x pa(y)]);                %other neighbors besides x and parent of y
                    D_prime(j) = 0;
                    if length(Ny)>0
                        for jj=1:length(Ny)
                            current_Ny = Ny(jj);
                            D_prime(j) = max(D_prime(j),D(current_Ny)+1);
                        end
                    end
                end
            end
        end
        Nx(remove_index)=[];
        D_prime(remove_index)=[];
        %% Choose y in Nx
        candidate_size = length(Nx);
        if candidate_size > 1
            horizontal_degree = zeros(1,candidate_size);
            for j=1:candidate_size
                y = Nx(j);
                Ny = tail(find(head==y));               %neighbor of y in G
                horizontal_degree(j)=length(find(U(Ny)==U(y)));
            end
            Nx_matrix = zeros(candidate_size,4);            %first col are the candidate nodes, last 3 cols are the constraint 
            Nx_matrix(:,1) = Nx';
            Nx_matrix(:,2) = U(Nx)';
            Nx_matrix(:,3) = D_prime';
            Nx_matrix(:,4) = horizontal_degree';
            Nx_matrix = sortrows(Nx_matrix,[-2 -3 -4]);        %in order of decreasing Uy, decreasing D'y, decreasing horizontal_degree
            y = Nx_matrix(1,1);
        else
            y = Nx;
        end
        %% Add edge into G_star 
        head_star = [head_star x y];
        tail_star = [tail_star y x];
        %% Modify T
        pa(x)=y;
        %% Update U values of x and descendant of x
        delta_x = U(y)+1-U(x);
        if delta_x ~= 0
            descendant = x;
            stack = x;
            while isempty(stack)
                new_nodes = find(pa==stack(1));
                descendant = [descendant new_nodes];
                stack = [stack new_nodes];
                stack(1)=[];
            end
            U(descendant)=U(descendant)+ delta_x;
        end
        %% update D values of y and ancestor of y
        delta_y = max(D(y),D(x)+1)-D(y);
        while delta_y ~= 0
            D(y)=D(y)+delta_y;
            delta_y = max(D(pa(y)),D(y)+1)-D(pa(y));
            y=pa(y);
        end
    end
end
T_star = sparse(head_star,tail_star,ones(length(head_star),1));

%% ---------------- Compute objective value fval ---------------- %%

nodes = unique(head_star);      %nodes in T_star       
node_number = length(nodes);    %number of nodes in T_star
%% Compute downward-distance D 
D = zeros(1,n);
visited = 0;
L = 0;
T_temp = T_star;
[T_temp_a, T_temp_b, w]=find(T_temp);
while visited<node_number-1
    one_vector_temp = ones(size(T_temp,1),1) ;
    degree = T_temp*one_vector_temp;
    leaf = find(degree==1);
    leaf = setdiff(leaf,intersect(leaf,root));
    D(leaf')=L;
    temp = ismember(T_temp_a,leaf);
    T_temp_a = T_temp_a(~temp);
    T_temp_b = T_temp_b(~temp);
    temp = ismember(T_temp_b,leaf);
    T_temp_a = T_temp_a(~temp);
    T_temp_b = T_temp_b(~temp);
    w = ones(length(T_temp_a),1);
    T_temp = sparse(T_temp_a, T_temp_b, w);
    visited = visited+length(leaf);
    L = L+1;
end
D(root)=L;

%% Compute degree of each node in T_star
degree = T_star*one_vector;

%% Compute objective value fval
fval = 2*D(root);
for i=1:node_number
    current_node = nodes(i);
    fval = fval+(degree(current_node)-2)*D(current_node);
end
                        
                    
        
    
    
    
    
    
