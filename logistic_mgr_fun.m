% assumes stress adjusted growth rates to follow a logistic function
% with a trade-off between Rmax and hT (Tolerance)
% 
% the trade-off curve is assumed to always go through the points [h0min  Rmin]
% and [hmax Rmax] so that the most tolerant species can barely survive at 0
% and the species with the highest growth rate always has the highest value
% of hT hTmax
% trade-off function is hT = a * R^n + b
%
% INPUT
%   h: elevation [double vector]
%   beta_r: Transition sharpness of r(h)
%   Rmax: values of Rmax of different species [double vector]
%   Rmaxmin: minimal value of Rmax of all species 
%   Rmaxmax: maximal value of Rmax of all species 
%   hTmax: maximal value of hT for trade-off 
%   out: location on domain at which most tolerant species can barely
%   survive (here: 0)
%   tradeoff: option variable [char] determining the trade-off
%
% OUTPUT
%   mgr: stress adjusted growth rates [matrix]: 
%        stress adjusted growth rates for every species at every elevation point
%   hT: hT values according to trade-off between Rmax and hT
%
% Mara Bannuscher, 15.04.2024, Matlab R2023a
% Torben Schucht, 13.01.2026

function [mgr, hT] = logistic_mgr_fun(h, beta_r, Rmax, Rmaxmin, Rmaxmax,hTmax, ...
    out, tradeoff, singleSp)

switch tradeoff

    case char('linear')
        n = 1;
    case char('upward')
        n = 1/4;
    case char('upwardstrong')
        n = 1/128;
    case char('downward')
        n = 2;
    case char('downwardstrong')
        n = 6;
    otherwise
        warning ('Please enter: linear, upward or downward')

end

hTmin = (log(Rmaxmin-1)/beta_r)+out;
    
% % trade-off parameters
a = (hTmin - hTmax)/(Rmaxmin^n-Rmaxmax^n);
b = hTmax - ((hTmin - hTmax)/(Rmaxmin^n-Rmaxmax^n))*Rmaxmax^n;

% TRADE-OFF
hT = a.*Rmax.^n+b;
    
if ~isempty(hT)
    %hT = hTmin;
    % stress adjusted growth rate FUNCTION
    mgr = Rmax./(1+exp(-beta_r.*(h-hT))); % logistic function 
end
        
end 

