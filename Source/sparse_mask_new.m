function [ sparse_out ] = sparse_mask_new( Img, Mask )
% Function is designed for calculating the weights on n-edges.
% Inputs
%     Img is the CT intensity
%     Mask is the region of interest
% Outputs
%     sparse_out records the weights of n-edges

Nsize = size( Img );

%% get the Adjacency Matrix, we consider 18-connected neighbors here.

% S_all = getAdjacencyMatrix(Nsize, 6);
S_all = getAdjacencyMatrix(Nsize, 18);
% S_all = getAdjacencyMatrix(Nsize, 26);
Ind_mask = find( Mask );
Len_ind = length(Ind_mask);
sparse_mask = S_all(Ind_mask, :);
sparse_mask = sparse_mask(:, Ind_mask);

Img_mask = Img( Mask );

[Ind_ii, Ind_jj, Dist] = find( sparse_mask);
clear sparse_mask;
weight = weight_fun(Img_mask(Ind_ii), Img_mask(Ind_jj), Dist);

sparse_out = sparse(Ind_ii, Ind_jj, weight, Len_ind, Len_ind);
sparse_out = sparse_out + sparse_out';

end



function weight = weight_fun(I1, I2, dist)
% get the edge (I1, I2) weight 
% maybe later on we can using other function April 3, 2015
%% early one
% dif = abs(I1-I2);
% weight =exp( -dif . distance );

% consider the distance of edge
dif = abs(I1-I2);
w2 = 1./(dist + 0.01);

if size(dif, 2) == size(w2, 2)
    weight = exp( -dif./w2 );
else
    weight = exp( -dif./w2' );
end

% if size(dif, 2) == size(w2, 2)
%     weight =w2.* exp( -dif );
% else
%     weight =w2'.* exp( -dif );
% end

end

