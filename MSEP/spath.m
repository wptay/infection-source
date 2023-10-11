function path= spath(A,a,b)
A=sparse(A);
path(1)=b;
A=double(A);
[d,pred] = shortest_paths(A,a);
i=b;
while pred(i) ~= 0
    path=[pred(i) path];
    i=pred(i);
end
end

