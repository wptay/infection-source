%SMALLW     Generate adjacency matrix for a small world network.
%
%   Input   n: dimension of matrix (number of nodes in graph).
%           k: number of nearest-neighbours to connect. Defaults to 1.
%           p: probability of adding a shortcut in a given row. Defaults to
%           0.1.
%
%   Output  A: n by n symmetric matrix with the attribute sparse.
%
%   Description:    Shortcuts are added to a kth nearest neighbour ring
%                   network with n nodes by calling the utility function
%                   short.m.
%
%   Reference:  D.J. Watts, S. H. Strogatz,
%               Collective Dynamics of Small World Networks,
%               Nature 393 (1998), pp. 440-442.
%
%   Example:    A = smallw(100,1,0.2);


function [explicit,infected]=function_infection_spread(G,M,MAX,source,explicit_percentage) 
[Ga,Gb,w] = find(G); % Ga and Gb represents linking relationship in the underlying Network

%Here starting the infection process
p = rand; %probability that a susceptible node get infected at the next time slot
%only the source is infected at t=0
infected=source; 
infected_count = length(source);

susceptible=[];
for j=1:length(Ga)
    if ismember(Ga(j),infected) && ~ismember(Gb(j),infected)
        susceptible(end+1) = Gb(j);
    end
end

while infected_count < M
    susceptible_temp = susceptible;
    for j=1:length(susceptible)
        if rand <= p 
            new_infected = susceptible(j);
            infected(end+1)=new_infected;
            susceptible_temp = setdiff(susceptible_temp, new_infected);
            new_susceptible=Gb(Ga==new_infected);
            new_susceptible =setdiff(new_susceptible, intersect(new_susceptible,infected));
            susceptible_temp = union(susceptible_temp,new_susceptible);
        end
    end
    susceptible = unique(susceptible_temp);
    infected_count = length(infected);
end
if infected_count > MAX
    infected = infected(1:MAX);
    infected_count = MAX;
end

if explicit_percentage == 0
    lower_limit = max(0,2-1/p);
    q=lower_limit+(1-lower_limit)*rand(1,infected_count);

    dice = rand(1,infected_count);        %throw dice to determine whether the infection of any node is observable 
    observable = ceil(q-dice);
    explicit = infected*diag(observable);
    explicit = setdiff(explicit,0);  
else
    explicit_count = ceil(infected_count*explicit_percentage);
    % use the function RANDPERM to create a random permutation of the integers 1 through infected_count, then pick the first explicit_count values 
    index = randperm(infected_count);
    explicit = infected(index(1:explicit_count));
    explicit = sort(explicit);
end