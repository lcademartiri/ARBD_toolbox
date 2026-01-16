function [degrees,gq,twobodycolls] = graphAnalysis(edges,qs)

	% --- building the adjacency matrix and the associated graph ---
	edges=unique(edges,'rows'); % isolate unique collisions
	edges=unique([edges;edges(:,2),edges(:,1)],'rows'); % add the mirrored collisions to ensure the resulting graph is undirected
	adjM=zeros(max(edges,[],'all'),max(edges,[],'all')); % initialize the adjacency matrix
	edgesind=sub2ind(size(adjM),edges(:,1),edges(:,2)); % find the linear indices corresponding to the pairs of colliders
	adjM(edgesind)=1; % populate the adjacency matrix
	G=graph(adjM); % build the graph
	% ---
	
	% --- eliminate isolated particles from consideration
	[bins,binsizes]=conncomp(G); % extract the graph size distribution
	idxgraph = binsizes(bins) >= 2; % identify subgraphs with more than 1 node
	G = subgraph(G, idxgraph); % eliminate 1-node graphs
	% ---

	% --- calculate the distribution of degrees of the nodes in the graph (i.e., coordination numbers) ---
	degrees=degree(G); % get all degrees
	tempedge=linspace(1,max(degrees)+1,max(degrees)+1)'; % define the bins
	[degrees,tempedge]=histcounts(degrees,tempedge); % get the distribution
	degrees=degrees'; % reformat the distribution to column vector
	% ---

	% --- calculate the characteristics (nodes#, edges#, and completeness) of all clusters ---
	[bins,~]=conncomp(G); % extract the subgraphs 
	for ig=1:max(bins) % loop over all clusters
		sG{ig,1}=subgraph(G,bins==ig); % build the cluster as a subgraph
		nodestemp=numnodes(sG{ig}); % count the nodes
		edgestemp=numedges(sG{ig}); % count the edges
		gq(ig,:)=[qs,nodestemp,edgestemp,double(edgestemp==(nodestemp*(nodestemp-1))/2)]; % compile step#, nodes#, edges#, and completeness
	end
	% ---
	
	% --- count and store two-body collisions as the total number of edges of the main graph 'G' ---
	twobodycolls=numedges(G);
	% ---
end
