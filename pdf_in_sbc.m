function pdf=pdf_in_sbc(POS,S,res)

	% radial pdf calculation in SBC conditions complete with peak position analysis
	% requires a particle positions in time [POS][N x 3 x T]
	% requires the size of the observation window S.br
	% requires the resolution of the pdf bins
	% the pdf analysis also yields the N(t) fluctuation spectrum for fuzzy tether systems
	% the pdf struct also contains the number density rho(t) fluctuation spectrum for fuzzy tether systems
	[Nall,~,T]=size(POS); 
	S.bv=(4/3)*pi*(S.br^3);
	
	fprintf('### Initializing: PDF Calculation ###\n');
	pdf.edges = [0;sort((2*S.br+eps:-res:0)')];
	pdf.edges2 = pdf.edges.^2;
	pdf.centers = pdf.edges(1:end-1,:)+0.5*diff(pdf.edges);
	pdf.vols = (4/3)*pi*(pdf.edges(2:end,1).^3-pdf.edges(1:end-1,1).^3);
	pdf.acc = zeros(size(pdf.centers,1),1);
	pdf.r_norm = pdf.centers / S.br; % Geometric Form Factor - Normalize distance by Box Radius (S.br)
	pdf.geom_factor = 1 - (3/4)*pdf.r_norm + (1/16)*pdf.r_norm.^3; % The Finite Volume Correction Polynomial
	pdf.geom_factor = max(0, pdf.geom_factor); % Clamp negative values to 0
	tStart=tic;
	counterstruct = struct('Stage','PDF Calculation');
	for it = 1:T
		p = squeeze(POS(:,:,it));
		pr = p(vecnorm(p,2,2)<S.br,:);
		pdf.N_t(it,1) = size(pr,1);
		pdf.rho_t(it,1) = pdf.N_t(it,1)/S.bv;
		sqdist = pdist(pr,'squaredeuclidean')';
		hc = histcounts(sqdist,pdf.edges2);
		denom = pdf.rho_t(it,1).*pdf.vols.*pdf.geom_factor;
		pdf.acc = pdf.acc + 2.*(hc(:)./denom)/pdf.N_t(it,1);
		progressUpdate(it, T, tStart, 100, counterstruct)
	end
	pdf.pdf = pdf.acc./T;
	
	
	% --- find peaks ---
	pdf.HRcenters = [0;sort((2*S.br+eps:-0.0005*S.rp:0)')];
	pdf.analysis = [pdf.HRcenters,interp1(pdf.centers,pdf.pdf,pdf.HRcenters,"pchip")];
	pdf.analysis(:,3) = pdf.analysis(:,2)-1;
	pdf.analysis(:,4) = pdf.analysis(:,3);
	pdf.analysis(:,5) = pdf.analysis(:,3);
	pdf.analysis(pdf.analysis(:,4)<0,4) = nan;
	pdf.analysis(pdf.analysis(:,5)>0,5) = nan;
	pdf.analysis(:,5) = abs(pdf.analysis(:,5));
	pdf.analysis(floor(end/2):end,:)=[];
	
	% find maxima & minima
	[pks, locs] = findpeaks(pdf.analysis(:,4));
	pdf.maxima=[locs,pks];
	pdf.maxima(pdf.maxima(:,2)<0.01*max(pdf.maxima(:,2)),:)=[];
	[pks, locs] = findpeaks(pdf.analysis(:,5));
	pdf.minima=[locs,pks];
	pdf.minima(pdf.minima(:,1)<=pdf.maxima(1,1),:)=[];
	pdf.minima(pdf.minima(:,2)<0.03*max(pdf.minima(:,2)),:)=[];
	pdf.maxima(:,1)=pdf.analysis(pdf.maxima(:,1),1);
	pdf.minima(:,1)=pdf.analysis(pdf.minima(:,1),1);
	fprintf('=== Completed: PDF Calculation ===\n');
end