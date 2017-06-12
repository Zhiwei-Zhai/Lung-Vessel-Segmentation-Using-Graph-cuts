addpath( genpath(pwd) );

T1=0.0009;
T2 = 0.8;
%% load data
Data = load('data.mat');
CT = im2ones(Data.CT);  % CT intensity
Strain = im2ones( Data.Strain );  %vesselness from strain energy filter: pxenhancement of https://github.com/ITKTools/ITKTools
AirwayDTF = Data.AirwayDTF;  

im = Strain;
mask1 = im > T1;
mask2 = im > T2;
mask = mask1 &( ~mask2);

%% Segment the lung vessels
tic, labels = Gc_Sparse(CT, Strain, mask, AirwayDTF); toc;

outImage = zeros( size( CT ) );     L_mask = ( labels > 0 ) | mask2;
outImage( L_mask ) = 1; 

%% show the results
sliceN = 10;

subplot(2,2,1), imshow(CT(:,:,sliceN)', []);
title('1: CT image');
subplot(2,2,2), imshow(Strain(:,:,sliceN)', []);
title('2: Vesselness');
subplot(2,2,3), imshow(AirwayDTF(:,:,sliceN)', []);
title('3: Airway DTF');
subplot(2,2,4), imshow(outImage(:,:,sliceN)', []);
title('4: Vessel Segmentation');

%% remove the path
rmpath( genpath(pwd) );