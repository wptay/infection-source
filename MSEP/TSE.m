function twoS_estimate = TSE(matrix, root)
one_source_logpermutation=[];
one_source_estimate=[];
global A;
A=matrix;
global N;
N = size(A) ;
N = N(1,1);
global t2_table;
global p2_table;
global p1_table;
global logp1_table;
global logp2_table;

global tt1_table;
global tt2_table;
global logpp1_table;
global logpp2_table;
global M2_table;
global qfunction_table;

t2_table = hashtable;
p2_table = hashtable;
t1_table = hashtable;
p1_table = hashtable;
logp1_table = hashtable;
logp2_table = hashtable;

tt1_table = hashtable;
tt2_table = hashtable;
logpp1_table = hashtable;
logpp2_table = hashtable;

logpermutation_table = hashtable;
M2_table = hashtable;
qfunction_table=hashtable;

% Functions for hashtable:
%     table = put(table, key, value);   - store a key/value pair
%     table = put(table, key);          - store just a key (to use as a Set)
% 
%     value = get(table, key);          - retrieve a stored value by key
% 
%     b     = has_key(table, key);      - test for key existence
% 
%     n     = count(table);             - get number of stored keys
% 
%     [k,v] = dump(table);              - get cell arrays of all keys and
%     values
%   
%             info(table);              - print table statistics and info

B=A;
C= ones(N,1) ;

D =B*C  ;   %get the degree of all nodes in B

while find(D(:,1))
    for i=1:N     
        if D(i)==1 && i~=root      % find the node i that has degree 1
            if has_key(t1_table, i);  
                t1_table =put(t1_table,i,1+get(t1_table, i));   %calculate the number of nodes in the subtree of node i, including itself
            else
                t1_table =put(t1_table,i,1);    %for the leaf nodes
            end
            if has_key(p1_table, i);
                p1_table =put(p1_table,i,get(t1_table, i)*get(p1_table, i));   %calculate the cummulative product of the subtree of node i, including itself
                logp1_table = put(logp1_table,i,log(get(t1_table,i))+get(logp1_table,i));
            else
                p1_table =put(p1_table,i,1);    %for the leaf nodes
                logp1_table = put(logp1_table,i,0);
            end
            j=find(B(i,:)==1);   %j would be the parent of node i
            t2_table =put(t2_table,[i j],get(t1_table, i));     % node i pass the number of nodes in its subree to its parent j
            p2_table =put(p2_table,[i j],get(p1_table, i));     % node i pass the cummulative product of its subtree to its parent j
            logp2_table = put(logp2_table,[i j],get(logp1_table, i));
            B(i,j)=0;       %get rid of the edge
            B(j,i)=0;
            if has_key(t1_table, j);  
                t1_table =put(t1_table,j,get(t1_table, i)+get(t1_table, j));    %node j sum up the number of nodes in all its children's subtrees
            else
                t1_table=put(t1_table,j,get(t1_table, i));
            end
            if has_key(p1_table, j);  
                p1_table =put(p1_table,j,get(p1_table, i)*get(p1_table, j));    %node j multiplies the cummulative product of all its children's subtrees
                logp1_table = put(logp1_table, j, get(logp1_table, i)+get(logp1_table,j));
            else
                p1_table=put(p1_table,j,get(p1_table, i));
                logp1_table=put(logp1_table,j,get(logp1_table,i));
            end
        end
    end
    D = B*C  ;
end

passdown(root,0);      % This is STEP 2, the root node passes the information down to the rest of the network. Here 0 indicates root node does not have a parent
%After passdown p1_table record the cummulative product of all the nodes'
%(excluding itself) subtree sizes considering itself as the source

logNf=0;
for i=1:N
    logNf=logNf+log(i);
end
logNminus2f=logNf-log(N)-log(N-1);

[key1,value1] = dump(logp1_table);

for i=1:length(value1)
    one_source_logpermutation=[one_source_logpermutation logNf-value1{i,1}-log(N)];
end

index_one_source_estimate=find(one_source_logpermutation==max(one_source_logpermutation));

if length(index_one_source_estimate) >= 2
    for i=1:length(index_one_source_estimate)
        one_source_estimate=[one_source_estimate key1{index_one_source_estimate(i),1}];
    end
else
    one_source_estimate=[one_source_estimate key1{index_one_source_estimate,1}];
    temp=sort(one_source_logpermutation,'descend');
    temp=temp(2);
    index_one_source_estimate=find(one_source_logpermutation==temp);
    one_source_estimate=[one_source_estimate key1{index_one_source_estimate(1),1}];
end

%STEP 3
D =A*C  ;   %get the degree of all nodes in A
B=sparse(A);
B=double(B);

leaf=find(D(:,1)==1);

%tt1_table [i j] recordes the subtree size rooted at node i directed away from node j
%tt2_table [i j k] recordes the subtree size rooted at node i directed away from the backbone
%logpp1_table [i j] recordes the cummulative product of the subtree size of all children of node i directed away from node j
%logpp2_table [i j k] recordes the cummulative product of the subtree size of all children of node i directed away from the backbone

for i=1:length(leaf)-1
    for j=i+1:length(leaf)
        path = spath(B, leaf(i), leaf(j));
        if ismember(one_source_estimate,path) 
            for k=1:length(path)
                switch k 
                    case 1
                        if ~has_key(tt1_table,[path(k) path(k+1)])
                            tt1_table =put(tt1_table,[path(k) path(k+1)],N-get(t2_table,[path(k+1) path(k)]));
                            logpp1_table =put(logpp1_table,[path(k) path(k+1)],get(logp1_table,path(k))-get(logp2_table,[path(k+1) path(k)]));
                        end
                    case length(path)
                        if ~has_key(tt1_table,[path(k) path(k-1)])
                            tt1_table =put(tt1_table,[path(k) path(k-1)],N-get(t2_table,[path(k-1) path(k)]));
                            logpp1_table =put(logpp1_table,[path(k) path(k-1)],get(logp1_table,path(k))-get(logp2_table,[path(k-1) path(k)]));
                        end
                    otherwise
                        if ~has_key(tt1_table,[path(k) path(k+1)])
                            tt1_table =put(tt1_table,[path(k) path(k+1)],N-get(t2_table,[path(k+1) path(k)]));
                            logpp1_table =put(logpp1_table,[path(k) path(k+1)],get(logp1_table,path(k))-get(logp2_table,[path(k+1) path(k)]));                   
                        end

                        if ~has_key(tt1_table,[path(k) path(k-1)])
                            tt1_table =put(tt1_table,[path(k) path(k-1)],N-get(t2_table,[path(k-1) path(k)]));
                            logpp1_table =put(logpp1_table,[path(k) path(k-1)],get(logp1_table,path(k))-get(logp2_table,[path(k-1) path(k)]));
                        end

                        if ~has_key(tt2_table, [path(k) path(k-1) path(k+1)])
                            tt2_table =put(tt2_table, [path(k) path(k-1) path(k+1)], N-get(t2_table,[path(k-1) path(k)])-get(t2_table,[path(k+1) path(k)]));
                            tt2_table =put(tt2_table, [path(k) path(k+1) path(k-1)], get(tt2_table,[path(k) path(k-1) path(k+1)]));
                            logpp2_table =put(logpp2_table, [path(k) path(k-1) path(k+1)], get(logp1_table,path(k))-get(logp2_table,[path(k+1) path(k)])-get(logp2_table,[path(k-1) path(k)]));
                            logpp2_table =put(logpp2_table, [path(k) path(k+1) path(k-1)], get(logpp2_table,[path(k) path(k-1) path(k+1)]));
                        end
                end
            end

            %% calculation of T_rho       
            T_matrix = zeros(max(path),max(path));
            for k=2:length(path)-2
                T_matrix(path(k),path(k+1))=get(tt2_table,[path(k) path(k-1) path(k+1)])+get(tt2_table,[path(k+1) path(k) path(k+2)]);
                if k+2 < length(path)
                    for x=k+2:length(path)-1
                        T_matrix(path(k),path(x))=T_matrix(path(k),path(x-1))+get(tt2_table,[path(x) path(x-1) path(x+1)]);
                    end
                end
            end
        
            %% prepare q value        
            q_matrix = zeros(max(path),max(path));
            for k=2:length(path)-1
                q_matrix(path(k),path(k))=1/get(tt2_table, [path(k) path(k-1) path(k+1)]);
            end

            for l=1:length(path)-3 %distance between v1 and vk
                for k=2:length(path)-1-l
                    q_matrix(path(k),path(k+l))=(q_matrix(path(k),path(k+l-1))+q_matrix(path(k+1),path(k+l)))/T_matrix(path(k),path(k+l));
                end
            end
        
            %% Calulate C value 
            for k=1:length(path)-1
                for x=k+1:length(path)
                    if ~has_key(logpermutation_table, [path(k) path(x)])
                        logpermutation = logNminus2f-get(logpp1_table,[path(k) path(k+1)])-get(logpp1_table,[path(x) path(x-1)]);
                        if x > k+1
                            logpermutation = logpermutation + log(q_matrix(path(k+1),path(x-1)));      
                            for y=k+1:1:x-1
                                logpermutation = logpermutation-get(logpp2_table,[path(y) path(y-1) path(y+1)]);
                            end
                        end
                        logpermutation_table=put(logpermutation_table,[path(k) path(x)], logpermutation);
                        logpermutation_table=put(logpermutation_table,[path(x) path(k)], logpermutation);
                    end
                end
            end
        end
    end
end

[key,value] = dump(logpermutation_table);

maximum=-inf;
for i=1:length(value)
    temp=value{i,1};
    if maximum < temp
        maximum=temp;
    end
end

twoS_estimate=[];
for i=1:length(key)
    if value{i,1}==maximum
        twoS_estimate=unique(union(twoS_estimate,key{i,1}));
    end
end

            

        
       
