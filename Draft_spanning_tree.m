%Find the minimum spanning tree and shows the graph
% m is the number of nodes
%method is the method utilized Prim or Kruskal
%Default parameters are m=10; method='prim';
function Draft_spanning_tree(m,method)

if nargin<1
    m=10;
    method='prim'; %'kruskal'
end

%Create Adjacency matrix
i=randi(m,1,m);
j=randi(m,1,m);
s=randi(10,1,m);
W=sparse(i,j,s,m,m);

%Undirected matrix
LW = tril(W + W');
LW = LW - diag(diag(LW));

view(biograph(LW,[],'ShowArrows','off','ShowWeights','on'))

%Finding the minimum spanning tree of graph   
[ST,pred] = graphminspantree(LW,m,'Method',method);

view(biograph(ST,[],'ShowArrows','off','ShowWeights','on'))
end