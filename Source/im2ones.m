function [ im ] = im2ones( imp )
% Function is used for rescale the imp to [0, 1];

N_min = min( imp(:) );
imT = imp - N_min;
imT = double(imT);
N_max = max( imT(:) );

im = double(imT./N_max);

end

