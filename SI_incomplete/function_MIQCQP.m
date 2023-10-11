function [fval,x]=function_MIQCQP(G, root)
alpha = 1;    %penalty for large D values (for the case x_ij = 1) 
%%----------------prepare needed parameter matrices                
[head,tail,w] = find(G);            %sparse representation, any edge if from tail to head
nodes = unique(head);               %array stores all nodes in increasing order
n = length(nodes);                  %number of nodes
m = length(head);                   %number of (directed) edges
nodes_index = zeros(1,max(nodes));  %store the index all nodes, for example, if nodes=[2 4 5], nodes_index=[0 1 0 2 3]
for i=1:n
    current_node = nodes(i);
    nodes_index(current_node) = find(nodes==current_node);
end
tail_matrix = zeros(m,n);   %extend the tail to a matrix
head_matrix = zeros(m,n);   %extend the head to a matrix
for i=1:n
    current_node = nodes(i);
    tail_matrix(tail==current_node,i)=1;
    head_matrix(head==current_node,i)=1;
end

%%----------------prepare input for optimization function 
dimention = m+n;
%------------------------------------H
H=zeros(dimention, dimention);
H(m+1:m+n,1:m)=tail_matrix';
H(1:m,m+1:m+n)=tail_matrix;
%------------------------------------f
f=zeros(dimention,1);
f(m+1:m+n,1)=-1*ones(n,1);       
f(m+nodes_index(root),1)=0;
f(m+1:m+n,1)=f(m+1:m+n,1)+alpha*ones(n,1);      %penalty for large D values
%------------------------------------Aeq
Aeq=zeros(n+1, dimention);
%constraint (1)
Aeq(1:n,1:m)=head_matrix';
%constraint (2)
Aeq(n+1,1:m)=ones(1,m); 
%------------------------------------beq
beq=zeros(n+1,1);          
%constraint (1)
beq(1:n)=ones(1,n);
beq(nodes_index(root))=0;
%constraint (2)
beq(n+1)=n-1;
%------------------------------------A
%constraint (6)
A=zeros(1,dimention);
A(1,1:m)=-tail_matrix(:,nodes_index(root))';
%------------------------------------b
%constraint (6)
b=-1;
%------------------------------------lb
lb=zeros(dimention,1);
%------------------------------------Q & l & qrl & qru %Quadratic Constraints (qrl <= x'Qx + l'x <= qru)
Q=cell(m+n,1);              %A cell array of double matrices, column orientated. Each cell is a constraint Q.
l=zeros(dimention,m+n);     %A double matrix, where each column is a constraint l vector.
qrl=-Inf(m+n,1);            %A double column vector, where each row is a constraint. 
qrl(m+1:m+n,1)=zeros(n,1);  %last n constraints serve as upper bound
qru=Inf(m+n,1);             %A double column vector, where each row is a constraint.
qru(1:m,1)=zeros(m,1);      %first m constraints serve as lower bound
for i=1:m
    %constraint (4)
    Q_temp = zeros(dimention,dimention);
    Q_temp(m+nodes_index(head(i)),i)=1;
    Q{i,1}=Q_temp;
    l(i,i)=1;
    l(m+nodes_index(tail(i)),i)=-1;
end
for i=1:n
    %constraint (5)
    current_node = nodes(i);
    Q_temp = zeros(dimention,dimention);
    index_neighbor = find(tail==current_node);
    for ii=1:length(index_neighbor)
        index_current_neighbor = index_neighbor(ii);
        current_neighbor = head(index_current_neighbor);  
        Q_temp(m+nodes_index(current_neighbor),index_current_neighbor)=1;
        l(index_current_neighbor,m+i)=1;
    end
    Q{m+i,1}=Q_temp;
    l(m+i,m+i)=-1;
end
%------------------------------------Integer Constraints
xtype = '';
for i = 1:m
    xtype = strcat(xtype,'I');
end
for i=1:n
    xtype = strcat(xtype,'C');
end
%%----------------Create OPTI Object
opts = optiset('maxnodes',100000);
Opt = opti('qp',H,f,'ineq',A,b,'eq',Aeq,beq,'lb',lb,'qcrow',Q,l,qrl,qru,'xtype',xtype,'options',opts);
%%----------------Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt);

