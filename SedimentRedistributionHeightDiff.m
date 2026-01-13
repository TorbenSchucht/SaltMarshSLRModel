% Sediment Redistribution based on height difference

%%% Input %%%

% hfull: Elevation profile
% D: Distance matrix between each point
% dx: Distance between two neighbor points
% lambdaFun: Sigmoidal Function to calculate lambda (width of kernel)
% lambda0: maximum width
% beta_l:   steepness of sigmoidal increase
% h_theta: elevation around which sigmoidal increase occurs

%%% Output %%%

% hfull: New Elevation profile

function hfull = SedimentRedistributionHeightDiff(hfull, D, dx, lambdaFun, lambda0, beta, h_thresh, lambdaMin)

initProfile = hfull;

% Height difference matrix
hdiff = hfull(:) - hfull(:)'; % hdiff_ij = "How much is i larger than j"
activation = hdiff>0; % Activation matrix

% Compute lambda for each point
lambda_vals = lambdaFun(lambda0, beta, h_thresh, hfull, lambdaMin); % Vector (N×1)
LambdaMatrix = lambda_vals'; % Lambda --> width of smoothing kernel

% Redistribution kernel for all pairs
K = (1./sqrt(2*pi*LambdaMatrix.^2)).*exp(-D.^2./(2*LambdaMatrix.^2)); % Gaussian kernel

% Apply redistribution
h_new = ((K .* hdiff)*dx).*activation;
h_cont = sum(h_new,2); % Shows how much each cell contributes to all cells, including itself (should be the old profile)
h_new = sum(h_new,1); % How much did each cell receive (from itself and all other cells; --> new profile)

hfull = h_new-h_cont'+initProfile; % Received - Contributed + InitialElevation

if abs(sum(hfull) - sum(initProfile)) > 1e-10 % Check that new and previous profile add to same sum
    disp( abs(sum(hfull) - sum(initProfile)))
end

clear D K % Remove for storage purposes
end