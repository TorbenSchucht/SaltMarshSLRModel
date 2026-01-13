% determines the depositional rate for the sediment accumulation

% Mara Bannuscher, 04.09.2024, Matlab R2023a
% Torben Schucht, 13.01.2026

% INPUT:
%   h: vector with elevation points corresponding to each grid point
%   cini: initial value for depositional rate, determining the steepness of
%         the function
%   d: initial value for combined losses
%
% OUTPUT:
%   c: depositial rate for each elevation point h
%       if c>0: deposition of sediment takes place
%       if c<0: erosion of sediment

function [c]  = deposition_erosion(h,co,d, beta_c)

c = co*log(1+exp(-beta_c * h))-d;

end


