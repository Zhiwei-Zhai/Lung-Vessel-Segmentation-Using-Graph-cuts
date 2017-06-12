function [ K ] = kernel_RBF_Mex_AirDTF( Img, Vsl, ADTF, Center, Sigma,  signal)
% the function was used for calcualting the weigt on t-edges. The intensity
% and vesselness of each voxel were considered. Airway distance transform
% was applied to compress the response of airway wall.
%   [ K ] = kernel_RBF_Mex( Img, Vsl, ADTF, Center, Sigma,  signal)
%   Inputs
%     Img is the CT ingensity
%     Vsl is the vesselness
%     ADTF is the Airway distance transform
%     Center is the Mu of img and Vsl
%     Sigma is the STD of img and Vsl
%   Outputs
%     K is the weight of t-edges

    w = 0.6; % the balance between 
    Prob = 0.95;
    P_const = sqrt( 2* log(2));
    Sigma = Sigma + 1e-6;

    % Prepare the energy for compressing airway wall.
    Mu_A = 1.2;  Sigma_A = 0.8;  w_A = 0.4;
    ADTF_cost = exp(-0.5*(ADTF-Mu_A).^2/Sigma_A.^2);
    ADTF_cost( ADTF==0 ) = 0; ADTF_cost( ADTF>3 ) = 0; 

    if signal ==1
        %% weights of t-edges between voxel nodes and sink node (background)
        beta = Center + P_const * Sigma;
        alpha = log(1/Prob - 1)./( P_const * Sigma );

        I_cost = Sigmoid_f(Img, alpha(1), beta(1) );
        Vsl_cost = Sigmoid_f(Vsl, alpha(2), beta(2));	

        K = w.*I_cost + (1-w ).*Vsl_cost + w_A.*ADTF_cost;
    else
        %% weights of t-edges between voxel nodes and source node (object)
        beta = Center - P_const * Sigma;
        alpha = - log(1/Prob - 1)./( P_const * Sigma );

        I_cost = Sigmoid_f(Img, alpha(1), beta(1) );
        Vsl_cost = Sigmoid_f(Vsl, alpha(2), beta(2));		

        K = w.*I_cost + (1-w ).*Vsl_cost - w_A.*ADTF_cost;
    end    

end



