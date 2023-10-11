function result = Choose_Sources_Geometric_Tree(matrix)
distance_between_two_sources=2;
one_side_distance=distance_between_two_sources/2;

A=matrix;
B=full(matrix);
N = size(B) ;
N = N(1,1);
C= ones(N,1) ;
D =A*C  ;

ref=(find(D==max(D)));
ref=ref(1);
% ref=ceil(N*rand(1,1));      %randomly choose one node and let's call it ref
[d dt pred] = bfs(A,ref);   

ref1=find(d==max(d));       %choose the node furthest away from ref, and let's call it ref1 (ref1 will be a leaf)
ref1=ref1(1);

[d dt pred] = bfs(A,ref1);  
ref2=find(d==max(d));       %choose the node furthest away from ref1, and let's call it ref2
ref2=ref2(1);
path= spath(A,ref1,ref2);   %considering the path from node ref1 to ref2
center_index=ceil(length(path)/2);
if center_index <= one_side_distance
    source=0;
else
    source=[path(center_index-one_side_distance) path(center_index+one_side_distance)];
end
result=source;