function path= spath(A,a,b)
%A=[0 1 0 0 0;1 0 1 1 0;0 1 0 0 0;0 1 0 0 1;0 0 0 1 0]; %test graph
A=sparse(A);
path(1)=b;
[d pred] = dijkstra_sp(A,a);
i=b;
while pred(i) ~= 0
    path=[pred(i) path];
    i=pred(i);
end
end

