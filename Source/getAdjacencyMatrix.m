function S_Matrix = getAdjacencyMatrix(Nsize, Pneigh)
% Function is for calculating the adjacency matrix for pixel/voxel in 2D/3D
% image. In 2D, the 4-connected and 8-connected types are
% available. In 3D, the 6-connected, 18-connected and 26-connected types
% are availabe.
% 
N = prod(Nsize);    % num of element in the array;
Ndim = numel( Nsize );

if Ndim == 3
    m = Nsize( 1 );    n = Nsize( 2 );     k = Nsize( 3 );
    
    switch lower( Pneigh )
        case {'6', 6}
            V1 = repmat( [ones(m-1, 1); 0], n*k, 1);    % for the left connect 1-2, 2-3; 
            V1 = V1(1: end-1);
            
            V2 = repmat( [ones(m*(n-1), 1); zeros(m, 1)], k, 1 );   % for below connect 1-5, 2-6, 3-7
            V2 = V2(1: end - m);
            
            V3 = ones( m*n*(k-1), 1);   % for under connect; 1-under1
            
            S_Matrix = sparse(1:(N-1), 2:N, V1, N, N)...
                + sparse(1:(N-m), (m+1):N, V2, N, N)...
                + sparse(1:(N - n*m), (m*n+1):N, V3, N, N);
        case {'18', 18}
            V1 = repmat( [ones(m-1, 1); 0], n*k, 1);  % for the left connect 1-2, 2-3;
            V1 = V1(1: end-1); % len(V1) = N-1
            S_Matrix = sparse( 1:(N-1), 2:N, V1, N, N );       clear V1;
            
            V2 = repmat( [ones(m*(n-1), 1); zeros(m, 1)], k, 1 );% for below connect 1-5, 2-6, 3-7
            V2 = V2(1: end - m); %length(V2) = N-m
            S_Matrix = S_Matrix + sparse( 1:(N-m), (m+1):N, V2, N, N );     clear V2
            
            V3 = ones( m*n*(k-1), 1); % length(V3) = m*n*(k-1)= N-m*n
            S_Matrix = S_Matrix + sparse( 1:(N-n*m), (m*n+1):N, V3, N, N ); clear V3;
            
            tmp_level = repmat( [ones(m-1, 1); 0], n-1, 1 );
            tmp_level = repmat( [tmp_level; zeros(m,1)], k, 1 ); 
            V5 = 1.4142*tmp_level(1: end-m);   clear tmp_level; %length(V5) = [m*(n-1)+m]*k -m = N-m
            V4 = V5(1: end-1);  % V4 = N-m-1
            S_Matrix = S_Matrix + sparse( 1:(N-m-1), (m+2):N, V4, N, N );   clear V4
            S_Matrix = S_Matrix + sparse( 2:(N-m+1), (m+1):N, V5, N, N);    clear V5;
            
            V7 = 1.4142*repmat( [ones(m-1, 1); 0], n*(k-1), 1); % length(V7) = m*n*(k-1)
            V6 = V7(1: end-1);  % length(V6) = m*n*(k-1) - 1
            S_Matrix = S_Matrix + sparse( 1:(N-m*n-1), (m*n+2):N, V6, N, N);    clear V6;
            S_Matrix = S_Matrix + sparse( 2:(N-m*n+1), (m*n+1):N, V7, N, N);    clear V7;
            
            V9 = 1.4142*repmat( [ones(m*(n-1), 1); zeros(m,1)], k-1, 1 );  % length of V9 is m*n*(k-1)
            V8 = V9(1: end-m);  %%length of V8 = m*n*(k-1) - m                         
            
            S_Matrix = S_Matrix + sparse( 1:(N-m*(n+1)), (1+m*(n+1)):N, V8, N, N);  clear V8;
            S_Matrix = S_Matrix + sparse( (m+1):(N-m*(n-1)), (1+m*n):N, V9, N, N);  clear V9;

        case {'26', 26}
            V1 = repmat( [ones(m-1, 1); 0], n*k, 1);  % for the left connect 1-2, 2-3;
            V1 = V1(1: end-1); % len(V1) = N-1
            S_Matrix = sparse( 1:(N-1), 2:N, V1, N, N );    clear V1;   
            
            V2 = repmat( [ones(m*(n-1), 1); zeros(m, 1)], k, 1 );% for below connect 1-5, 2-6, 3-7
            V2 = V2(1: end - m); %length(V2) = N-m
            S_Matrix = S_Matrix + sparse( 1:(N-m), (m+1):N, V2, N, N ); clear V2;
            
            V3 = ones( m*n*(k-1), 1); % length(V3) = m*n*(k-1)= N-m*n
            S_Matrix = S_Matrix + sparse( 1:(N-n*m), (m*n+1):N, V3, N, N ); clear V3;

            tmp_level = repmat( [ones(m-1, 1); 0], n-1, 1 );
            tmp_level = repmat( [tmp_level; zeros(m,1)], k, 1 ); %length(tmp_level) = m*n*k = N
            V5 = 1.4142*tmp_level(1: end-m);   %length(V5) = [m*(n-1)+m]*k -m = N-m
            V4 = V5(1: end-1);  % V4 = N-m-1
            S_Matrix = S_Matrix + sparse( 1:(N-m-1), (m+2):N, V4, N, N );   clear V4;
            S_Matrix = S_Matrix + sparse( 2:(N-m+1), (m+1):N, V5, N, N);    clear V5;
            
            V7 = 1.4142*repmat( [ones(m-1, 1); 0], n*(k-1), 1); % length(V7) = m*n*(k-1)
            V6 = V7(1: end-1);  % length(V6) = m*n*(k-1) - 1
            S_Matrix = S_Matrix + sparse( 1:(N-m*n-1), (m*n+2):N, V6, N, N);    clear V6;
            S_Matrix = S_Matrix + sparse( 2:(N-m*n+1), (m*n+1):N, V7, N, N);    clear V7;
            
            V9 = 1.4142*repmat( [ones(m*(n-1), 1); zeros(m,1)], k-1, 1 );  % length of V9 is m*n*(k-1)
            V8 = V9(1: end-m);  %%length of V8 = m*n*(k-1) - m
            S_Matrix = S_Matrix + sparse( 1:(N-m*(n+1)), (1+m*(n+1)):N, V8, N, N);  clear V8;
            S_Matrix = S_Matrix + sparse( (m+1):(N-m*(n-1)), (1+m*n):N, V9, N, N);  clear V9;
            
            V11 = 1.7321*tmp_level(1: end-m*n);  %clear tmp_level; %length(V11) = N-m*n = m*n*(k-1)
            V10 = V11(1: end-m-1); %length(V10) = len(V11)-m-1 = N-m*n-m-1
            S_Matrix = S_Matrix + sparse( 1:(N-m*(n+1)-1), m*(n+1)+2:N, V10, N, N); clear V10;
            S_Matrix = S_Matrix + sparse( m+2:N-m*(n-1)+1, m*n+1:N, V11, N, N); clear V11;
            
            V12 = 1.7321*tmp_level(1: end-m*(n+1));
            V13 = tmp_level(1:end-m*n-1); clear tmp_level;      
            
            S_Matrix = S_Matrix + sparse( 2:N-m*(n+1)+1, m*(n+1)+1:N, V12, N, N);   clear V12;
            S_Matrix = S_Matrix + sparse( m+1:N-m*(n-1)-1, m*n+2:N, V13, N, N); clear V13;
%             clear V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13;
        otherwise
            error('Sparse adjance matrix: Neighbor_type error', 'Unknown neighbor p for 3D (should be either 6 18 26)'); %#ok<CTPCT>
    end

elseif Ndim == 2
    m = Nsize( 1 );    n = Nsize( 2 );
    
    switch lower( Pneigh )
        case {'4', 4}
            % % 1-off diagonal elements
            V1 = repmat([ones(m-1,1); 0],n, 1);
            V1 = V1(1:end-1); % remove last zero
            % % n-off diagonal elements
            V2 = ones(m*(n-1), 1);
            
            % % get the upper triangular part of the matrix     %[1  2  3  4 ]
            S_Matrix = sparse(1:(N-1), 2:N, V1, N, N)...      %[5  6  7  8 ]
                + sparse(1:(N-m),(m+1):N, V2, N, N);             %[9  10 11 12]
        case {'8', 8}                                           %[13 14 15 16]
            V1 =  repmat([ones(m-1,1); 0],n, 1) ; % for the left connect 1-2, 2-3, 3-4; 
            V1 = V1(1:end-1);
            
            V2 =  ones(m*(n-1), 1) ;   % for below connect 1-5, 2-6, 3-7 4-8
            
            V4 = 1.4142*repmat([ones(m-1,1);0], n-1, 1) ;% for the below left connect, 1-6, 2-7, 3-8
            V3 = 1.4142*V4(1: end-1);      % the same with V1 ??????

            S_Matrix = sparse(1:(N-1), 2:N, V1, N, N)... 
                + sparse(1:(N-m), (m+1):N, V2, N, N)...
                + sparse(1:(N-m-1), (m+2):N, V3, N, N)...
                + sparse(2:(N-m+1), (m+1):N, V4, N, N);
            
        otherwise 
            error('Sparse adjance matrix: Neighbor_type error','Unkown neighbor p for 2D, (should be either 4 or 8)')
    end
else
    error('Dimension is not correct:', 'the dimension should be 2 or 3');
end

end
