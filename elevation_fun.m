% different versions for initial elevation
%
% INPUT
%   L: length of domain [float]
%   np: number of grid points [int]
%   dx: length of subintervals [float]
%   h_fun: decision variable, determining the initial elevation [char]
%   borders: variable, determining the border effects ('reflecting or
%   window' [char]
%                                                                       
%
% OUTPUT
%   h: initial elevation [vector]
%
% Plateau with increase + Reflecting boundaries is used in Ch.4
% Other options are experimental
%
% Mara Bannuscher, 09.05.2024, Matlab R2023a
% Torben Schucht, 13.01.2026


function [h] = elevation_fun(L,np,dx, h_fun, borders)

switch borders 

    case {'hanning', 'zero-padding'}
        
        switch h_fun
      
            case 'linear'
                h = linspace(-L, L-dx, 2*np);

            case 'less steep'
                h = linspace(-L/4, L/4 -dx, 2*np);
            
            case 'one tidal creek'
                tx = 0.5;    % Position of tidal creek
                tl = 0.3;    % length of tidal creek
                th = -0.2;  % depth of tidal creek (height of bottom)
                xx = linspace(-L, L-dx, 2*np);
                tpos = xx>=tx-1/2*tl & xx<=tx+1/2*tl; % get indices of xx
                h = linspace(-L, L-dx, 2*np);
                h(tpos) = th;  

            case 'island'

                % length of the plateu in percent area of the domain
                plateau = 0.2; 

                % double domain for window
                xx = linspace(-L,L-dx, 2*np);  

                % find first index of xx where xx >=-plateau/2
                plateau_first_index = find(xx>=-plateau/2,1,'first');
                plateau_last_index = find(xx>=plateau/2,1,'first');
                h1 = 2*linspace(-L,-plateau/2,plateau_first_index)+1.5;
                h2 = ones(1,plateau_last_index-plateau_first_index);
                h3 = -2*linspace(plateau/2,L,2*np-plateau_last_index)+1.5;
                h = [h1,h2,h3];


        end 

    case 'reflecting'

        switch h_fun

            case 'linear'
                h = linspace(-1.5*L,1.5*L-dx,3*np);
                
            case 'plateauWithIncrease'
                plateau = -1*ones(1,0.5*np);
                increase = linspace(-1,0.5*L-dx , 0.5*np);
                h = [plateau, increase];
                h = [fliplr(h),h,fliplr(h)];
             
            case 'less steep'
                h = linspace(-L,L-dx,3*np);
            
            case 'one tidal creek'
                tx = 0.5;    % Position of tidal creek
                tl = 0.3;    % length of tidal creek
                th = - 0.2;  % depth of tidal creek (height of bottom)
                xx = linspace(-1.5*L, 1.5*L-dx, 3*np);
                tpos = xx>=tx-1/2*tl & xx<=tx+1/2*tl; % get indices of xx
                h = xx;
                h(tpos) = th; 


            case 'island'
                % length of the plateu in percent area of the domain
                plateau = 0.5; 

                % triple domain for reflection
                xx = linspace(-1.5*L, 1.5*L-dx, 3*np);  

                % find first index of xx where xx >=-plateau/2
                plateau_first_index = find(xx>=-plateau/2,1,'first');
                plateau_last_index = find(xx>=plateau/2,1,'first');
                h1 = 2*linspace(-1.5*L,-plateau/2,plateau_first_index)+1.5;
                h2 = ones(1,plateau_last_index-plateau_first_index);
                h3 = -2*linspace(plateau/2,1.5*L,3*np-plateau_last_index)+1.5;
                h = [h1,h2,h3];

                
                case 'emergingIsland'
                % length of the plateu in percent area of the domain
                plateau = 0.2;
                
                % Elevation of plateau
                maxElevation = 0.3;

                % triple domain for reflection
                xx = linspace(-1.5*L, 1.5*L-dx, 3*np);  

                % find first index of xx where xx >=-plateau/2
                plateau_first_index = find(xx>=-plateau/2,1,'first');
                plateau_last_index = find(xx>=plateau/2,1,'first');
               
                h1 = linspace(-5 ,maxElevation , plateau_first_index);
                h2 = ones(1,plateau_last_index-plateau_first_index).*maxElevation;
                h3 = linspace(maxElevation, -5, 3*np-plateau_last_index);
                h = [h1,h2,h3];

        end 
end

