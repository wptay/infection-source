function passdown(node_name,parent)
global t2_table;
global p2_table;
global p1_table;
global logp1_table;
global logp2_table;
global A;
global N;
a=node_name;
C= ones(N,1) ;
D =A*C  ;
if parent==0       %check whether this node is root node
    child=find(A(a,:));     %find all the children of node a
else
    child=setdiff(find(A(a,:)),parent);      %find all the children of node a (all the neighbours of node a except for its parent)
end
for i=1:length(child)   %every loop passes the information to one of its children
    j=child(i);  
    t2_table =put(t2_table,[a j],N-get(t2_table, [j a]));
    p2_table =put(p2_table,[a j],get(p1_table, a)*get(t2_table, [a j])/get(p2_table, [j a]));
    logp2_table =put(logp2_table,[a j],get(logp1_table, a)+log(get(t2_table, [a j]))-get(logp2_table, [j a]));
    p1_table =put(p1_table,j,get(p2_table, [j a])*get(p2_table, [a j])/get(t2_table, [j a]));
    logp1_table =put(logp1_table,j,get(logp2_table, [j a])+get(logp2_table, [a j])-log(get(t2_table, [j a])));
    if D(j)~=1      %stop pass the information further if node j is a leaf
        passdown(j,a)
    end
end
