function [ Labels ] = Gc_Sparse(varargin )
% Do graph cut only in the region of interest, Construct the sparse graph
%  only in the mask region of the whole image. 
%   Detailed explanation goes here
%   1. [ Labels ] = KernelGc_Sparse_Shape_Journal( CT, Strain, mask )
%   2. [ Labels ] = KernelGc_Sparse_Shape_Journal( CT, Strain, mask, AirwayDTF )
%   Input:
%       CT
%		Strain
%		mask
%		AirwayDTF
%   Output:
%       Labels, the label result of the 
%       


Ninput = numel( varargin );
if ~isequal( Ninput, 3) && ~isequal(Ninput, 4)
    error('the input is error');
end
 
Fim = varargin{ 1 };
Vsl = varargin{ 2 };
Mask = varargin{ 3 };
if isequal( Ninput, 4)
	AirwayDTF = varargin{ 4 };
	RegionA = AirwayDTF( Mask );
else
    AirwayDTF = zeros(size(Fim));
	RegionA = AirwayDTF( Mask );
end


if ~isequal(size(Fim), size(Mask))
    error('the size of image and size of mask are not the same');
end

Region = Fim( Mask );
RegionV = Vsl( Mask );      %clear Vsl;
Nsize = size( Fim );


%% means as threshold
Thr1 = 1.04*mean(Region(:));     
Mu1 = [mean(Region(Region < Thr1)), 1.0*mean( Region(Region >= Thr1)) ];    
Sigma1 = [std(Region(Region < Thr1)), std(Region(Region >= Thr1)) ];

Thr2 = mean(RegionV(:)); 
Mu2 = [mean(RegionV(RegionV < Thr2)), mean(RegionV(RegionV >= Thr2))];
Sigma2 = [std(RegionV(RegionV < Thr2)), std(RegionV(RegionV >= Thr2))];

%% the paramenter important to set
k = 2; %k cluester to be classified 
gamma = 0.01;


Dc = zeros(k, length(Region));


for ci = 1 : k % initial the Dc
    % calculate the weight on t-edges (connect between source/sink node and voxel nodes)
    K = kernel_RBF_Mex_AirDTF(Region, RegionV, RegionA, [Mu1(ci), Mu2(ci)], [Sigma1(ci), Sigma2(ci)], ci);       
    Dc(ci,:)=1 - K;
end
clear K;
Sc = gamma*(ones( k ) - eye( k ));

% calculating the N-edges (connection between neighboring voxels)
sparSc = sparse_mask_new(Fim, Mask); 
clear Fim Vsl AirwayDTF;
%% graph cut module 
gch = GraphCut('open', Dc, Sc, sparSc);
[gch, L] = GraphCut('swap',gch);
[gch] = GraphCut('energy', gch);
gch = GraphCut('close', gch);

%% output
Labels = zeros( Nsize );
Labels(Mask) = L;

end