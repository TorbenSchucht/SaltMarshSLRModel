% Trait-based stress gradient model with sedimentation and redistribution

cd("C:\Users\torbe\Nextcloud\Shared\PhD Project\GitCloud\Torben\Second Paper\Code\")

% needed functions:
%   - logistic_mgr_fun
%   - elevation_fun
%   - deposition_erosion
%   - SedimentRedistributionHeightDiff
%   - plotVerticalChange (optional)

% Mara Bannuscher, 04.09.2024, Matlab R2023a
% Torben Schucht, 13.01.2026

clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS AND MODEL SPECIFICATIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General simulation properties (space,time)
Tsteps = 50000; % Max timesteps

% initial elevation: linear, one tidal creek, island, emergingIsland, less steep, 'plateauWithIncrease'
h_fun = 'plateauWithIncrease';

% Name of simulation
simname = 'xx';

L = 2;             % length of domain
n = 2^9;           % Spatial resolution / number of grid points
halfExt = 'on';     % Only use half as reflective extension --> Save computation time
%% Plotting of results
visual = true; % Make video
visualsteps = 10; % Timesteps after which a frame is generated
communityPoints = [100, 300, 500]; % Timepoints when community is plotted
tthres = 100;      % start of plot of Effective species richness


%% Community properties

% trade-off: linear, upward, upwardstrong, downward, downwardstrong
tradeoff = 'upwardstrong';

m = 20;            % number of species
singleSp = 'n';    % only generate single species
inidens = 0.01;    % initial density of every species at every point
sigma = 0.02;      % Dispersal width
thres = 0.01;      % threshold density under which species become extinct
out = 0;           % point of the domain at which the most tolerant species can barely survive

% parameters for stress-adjusted growth rates
Rmaxmax = 3;    % maximal maximal growth rate
Rmaxmin = 1.5;  % minimal maximum growth rate
hTmax = 0.7;    % highest inflection point of growth curves r(h)
beta_r = 15;    % Transition sharpness r(h)

%% Erosion parameters
d = 0.01;    % Combined losses (autocompaction and decomposition)


%% Sedimentation parameters
c0 = 0.1;          % Sedimentation intensity
beta_c = 2.5;      % Steepness sedimentation function

%% Sediment redistribution
lambda0 = 0.05;             % Redistribution intensity
beta_l = 2;                 % beta_lambda; Transition sharpness of redistribution function
h_theta = -0.5;             % Inflection height of redistribution function
lambdaMin = 0.01;           % Minimum redistribution

%% Sea-Level-Rise parameters

% a sea level rise in the model is equal to lowering the complete elevation
% regime by the factor sea_level_rise in each generation
sea_level = 'rising';   % rising or constant
sea_level_rise = 0.002;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION OF MODEL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize space

% Initialize spatial grid
dx = L/n;         % length of subintervals
x = linspace(-L/2,L/2-dx,n);    % domain/viable habitat
% triple domain for reflection
xx = linspace(-1.5*L, 1.5*L-dx, 3*n);
[XI, XJ] = meshgrid(xx, xx); % Create distance matrices
D = abs(XI - XJ);
% indices of the domain
domain = find(xx>=-L/2,1,'first'):1:(find(xx>=L/2,1,'first')-1);
if strcmp(halfExt, 'on') % If Extension is set to half size
    ext = n/2;
    Dcompressed = D(min(domain)-ext:max(domain)+ext,min(domain)-ext:max(domain)+ext);
end

% initialize elevation profile
h = elevation_fun(L,n,dx,h_fun,'reflecting');
hplot = h;   % safe initial value for end plot
h_star = -log(d/c0-1)/beta_c; % approx. Equilibrium elevation, which marsh will grow towards

%% Species initialization
% Generate Rmax values
Rmax = linspace(Rmaxmin,Rmaxmax,m)';  % Rmax values for all species
Rmaxplot = Rmax;                      % safe initial values for plot
% safe initial mgr and hT for all species for plot
[mgrplot,hT] = logistic_mgr_fun(xx, beta_r, Rmax, Rmaxmin, Rmaxmax, ...
    hTmax, out, tradeoff);

hTplot = hT;  % safe initial growth functions

if singleSp == 'y'
    mgrplot = max(mgrplot);
    hT = hT(1);
    m = 1;
end


% Dispersal kernel
LAPLACE = (1/sqrt(2*sigma^2)).* exp(-sqrt(2/(sigma^2)) .* abs(xx));
FLAPLACE = fft(LAPLACE); % Fourier transform of Laplace kernel

% initial population density Ninit for every species is iniden
% at every grid point for every species
Ninit = repelem(inidens,length(xx));

% create matrix to store densities for every species
% each row represents one species
N = repmat (Ninit,m,1);

N_id = (1:m)';  % vector with species ID

%% Initialize sedimentation / redistribution / SLR

% determine depositional rate for every start elevation
c = deposition_erosion(h, c0, d, beta_c);

% Redistribution
redistributionKernel = @(d, lambda) exp(-abs(d) ./ lambda) ./ (2 * lambda);

% Sigmoidal increase
lambdaFun = @(lambda0, beta_l, h_theta, hlocal, lambdaMin) lambda0 .* exp(beta_l * (hlocal - h_theta)) ./ (1 + exp(beta_l * (hlocal - h_theta)))+lambdaMin;

% SLR
sea_level_rise_vector = repmat(sea_level_rise,1,n);

%% Plot initialization

% Save community points for plot
hProfilePoints = [1,communityPoints(1),communityPoints(2),communityPoints(3)];

% prepare plots
if visual == true
    figure (1)
    tiledlayout(2,1)
end

% determine colourmap for density plots
cmap = colormap(jet(m));
cmapplot = cmap;   % safe original values for mgr plot

% default plot properties
set(0, 'DefaultAxesFontName', 'Times New Roman');  % For axes labels and titles
set(0, 'DefaultTextFontName', 'Times New Roman');  % For text objects
set(groot, 'defaultColorbarFontName', 'Times New Roman') % for colorbars
set(groot, 'DefaultLegendFontName', 'Times New Roman') % for legends
set(groot, 'DefaultColorbarFontSize', 15);
set(groot, 'DefaultLegendFontSize', 15);
set(groot,'DefaultTextFontSize',15)

% if the initial elevation isn't linear: create c values for a linear
% elevation for the end plot
if ~strcmp(h_fun,'linear')
    cplot = deposition_erosion(xx, c0, d, beta_c);
else
    % safe initial values for plot
    cplot = c;
end

%% Predefine some matrices for result saving

dens_max_iter = NaN(m,Tsteps+1);  % matrix for iteration
dens_max = NaN(m,Tsteps+1);       % end matrix (also contains extinct species)
dens_max_iter_ele = NaN(m,Tsteps+1);
dens_max_ele = NaN(m,Tsteps+1);   % determine the elevation at which each...
                                  % ...species reaches their maximum
dens_max_iter(:,1) = repmat(inidens,m,1); % store beginning densities:

%% Lower boundary of the salt marsh

xx_domain = xx(domain); % only evaluate domain (without borders)
% the index therefore refers to the index inside of the domain

front = NaN(1,Tsteps); % initialize vector for salt marsh front
front_index = find(h(domain)>0,1); % find index of front
front(1) = xx_domain(front_index); % initialize first value

% initialize vector for speed of front
velocities = NaN(1,Tsteps);

% only evaluate h inside of the domain
h_domain = h(domain);


%% Vertical accretion
va = NaN(Tsteps, n); % initialize matrix

%% Redistribution contribution
         redistributionChange = NaN(Tsteps, n); % Tracks changes in redistribution

         % Track elevation change and vertical change at front
         totalelevationChange = NaN(Tsteps,n);
         elevationChangeFront = NaN(1,n);

%% create video object

if visual == true
    videofilename = simname;
    v = VideoWriter(videofilename, 'MPEG-4');
    v.FrameRate = 10;   % picture rate: 10 per second
    open(v)
end

%% Iteration for Tsteps generations

% initialize matrix for effective species richness
SH = NaN(Tsteps, n);

% in every generation (year) we first assume growth and dispersal to
% happen, then calculate the sedimentat accumulation of this generation and
% apply the elevation change to the next generation

% Predefine biomass
biomassMeasureInterval = 20; % How often biomass should be measured across the gradient
biomass = zeros(ceil(20/Tsteps), length(h));
bIdx = 0; % Index for biomass saving
subplotIdx = 0;
% Predefine elevation profiles (for travelling waves)
htime = NaN(length(hProfilePoints),length(h));
htime(1,:) = h; % Save initial profile
hidx = 1;

for t=1:Tsteps % iterate through all generations/year
    
    % initialize competition vector before each growth phase
    cv = sum(N,1);
    
    % get current stress-adjusted growth rate and hT of all surviving species
    [mgr,hT] = logistic_mgr_fun(h, beta_r, Rmax, Rmaxmin, Rmaxmax, hTmax, out, tradeoff);
    if singleSp == 'y'
        mgr = max(mgr);
        hT = hT(1);
        m = 1;
    end
    %% if total species density is below threshold density:
    
    index_extinct = sum(N(:,domain),2)<thres; % index of extinct species
    
    % get species ID from extinct species
    id_extinct = N_id(index_extinct);
    
    % delete some properties/traits
    mgr(index_extinct,:) = [];    % delete mgr from species
    Rmax(index_extinct) = [];     % delete Rmax value from species
    hT(index_extinct) = [];       % delete hT value from species
    cmap(index_extinct,:) = [];   % delete color of species

    
    % save maximum (or sum) species densities over time
    dens_max(id_extinct,:) = dens_max_iter(index_extinct,:);
    dens_max_ele(id_extinct,:) = dens_max_iter_ele(index_extinct,:);
    % delete species in dens_max_iter matrix
    dens_max_iter(index_extinct,:) = [];
    dens_max_iter_ele(index_extinct,:) = [];
    
    % delete IDs from extinct species
    N_id(index_extinct) = [];
    
    % species is extinct -> delete species
    N(index_extinct,:) = [];
    
    
    %% calculate new population density
    
    % calculate new population density along the domain:
    Ni = mgr.*N.*exp(-cv); % growth function
    
    % densities in domain after growth:  
    % with reflecting boundaries:
    Ni_domain = Ni(:,domain);
    Ni = [fliplr(Ni_domain),Ni_domain,fliplr(Ni_domain)];
    
    % Fourier transform of new population density:
    FN = fft(Ni,[],2);  % along second dimension -> per row
    
    % inverse Fourier transform:
    N = dx*real( fftshift( ifft( FN.*FLAPLACE, [], 2 ),2 ) );
    
    
    %%  Sedimentation
    
    % linear relationship of Sedimentation S and biomass
    % biomass is approxied by total species density cv
    S = c.*cv;      % c is elevation-dependend
    
    % safe vertical accretion in matrix
    va(t,:) = S(domain);
    
    % determine new elevation regime
    h_domain = h_domain + S(domain);          % Sediments are added to elevation
    h = [fliplr(h_domain), h_domain, fliplr(h_domain)];
    
    % reflecting boundary condition:
    h_domain = h(:,domain);
    h_b4smooth = h_domain;
    
    %% Sediment redistribution
    
    if strcmp(halfExt, 'on') % Half extension
        hcompressed = h(min(domain)-ext:max(domain)+ext);
        
        % Apply redistribution to compressed space
        [hcompressed] = SedimentRedistributionHeightDiff(hcompressed, Dcompressed, dx, lambdaFun, lambda0, beta_l, h_theta, lambdaMin);        
        hcompressed = hcompressed(1+ext:length(hcompressed)-ext);
        h = [fliplr(hcompressed), hcompressed, fliplr(hcompressed)];
    else
        % Apply redistribution
        h = SedimentRedistributionHeightDiff(h, D, dx, lambdaFun, lambda0, beta_l, h_theta, lambdaMin);        
    end
    h_domain = h(domain);
    h = [fliplr(h_domain), h_domain, fliplr(h_domain)]; % Reflective boundaries
    
    h_aftersmooth = h_domain;
    
    % Optional plot to plot contribution of redistribution, sedimentation
    % and SLR:
    %plotVerticalChange(redistributionChange, va, totalelevationChange, front_index, t,x)
    
    %% Determine deposition for next year
    
    % new values for c:
    % c is based on the elevation regime of the current generation and
    % is used to determine the depositional rate for the next
    % generation
    
    c = deposition_erosion(h,c0,d, beta_c);
    
    % if density drops below zero at one point,
    % manually set it to 0 (for model performance)
    N(N<0) = 0;
    
    % check if whole gradient is submerged, or all species died
    % --> collapse
    if max(h) < 0 || size(N,1) == 0
        break
    end
    
     % get maximum density and index of maximum density:
    [densmax,idxmax] = max(N(:,domain),[],2);
    % maximum density of each species in domain
    dens_max_iter(:,t+1) = densmax;
    % get elevation of each species maximum elevation
    dens_max_iter_ele(:,t) = h_domain(idxmax);
    
    %% salt marsh front
    
    % front is defined as the point, where elevation = 0, which is the
    % elevation at which the most tolerant species can barely survive
    
    xx_domain = xx(domain); % only evaluate domain (without borders)
    
    % when looking at an island, we determine both fronts
    
    front_index = find(h(domain)>0,1,'first'); % find index of front
    if ~isempty(front_index)
        front(t+1) = xx_domain(front_index);
    end
    
    % determine velocity
    if t >1
        velocities(t-1) = front(t) - front(t-1); % horizontal velocity
    end
    
    %% Determine mean and max elevation
    
    % only evaluate h inside of the domain
    h_domain = h(domain);
    
    if isempty(front_index) || front_index == 1 % Front reached the border
        break
    end
    if isempty(N) % No species left
        break
    end
    
    % apply sea level change for next generation
    if strcmp(sea_level,'rising')
        h = h - repmat(sea_level_rise_vector,1,3);
        h_domain = h(domain);
    end
    
    %% Calculate effective species richness
    if t>tthres
        
        % initialize all species in domain
        Nd = N(:,domain);
        p = Nd./sum(Nd,1);           % relative abundance matrix
        
        % calculate Effective species richness:
        SH(t ,:) = exp(-sum(p.*log(p+(10^-13)))); % Add fraction of 10^-13 to avoid log(0)
        SH(t,sum(Nd,1)<1e-14) = 1; % Set richness to 1, if total biomass is extremely low at that position --> avoid numerical artefacts
    end
    
        %% Density and elevation plot
    
    % if visual is true: plot at every given generation vissteps:
    if visual == true && rem(t, visualsteps) == 0
        
        % plot densities
        figure(1),set(gcf, 'Name', 'Iteration', 'NumberTitle', 'off');
        
        % write video
        if visual == true
            writeVideo(v,getframe(gcf));
        end
        
        % plot current densities of all non-extinct species
        nexttile(1)
        pl = plot(xx,N,LineWidth=1.5);
        cbar = colorbar('Direction','reverse');
        cbar.Colormap = cmapplot; % use original values of cmapplot
        cbar.YTick = [0.05, 0.5, 0.95];
        cbar.TickLabels = {'high', 'medium', 'low'};
        cbar.Title.String = "Tolerance";
        cbar.Title.FontName = 'Times New Roman';
        cbar.FontSize=14;
        set(pl, {'color'}, num2cell(cmap,2)) % always use same colours
        xlim ([-1 1])         % set x-axis limits to domain
        
        ylabel('Density $N_{i}$', Interpreter='latex', FontSize=14)
        
        xticks([-1 -0.5 0 0.5 1])
        xticklabels({'-1','-0.5','0','0.5','1'})
        yticks([0 0.5 1 1.5])
        yticklabels({'0','0.5','1', '1.5'})
        ylim([0 1.5])         % set y lim so that plot size doesn't change
        
        title('Species densities', FontSize=15)
        set(gca, "FontSize",14)
        subtitle(['after ', num2str(t),' generations'], FontSize=14)
        drawnow
        
        % plot elevation change
        nexttile(2)
        plot(xx,h,LineWidth=2)
        xline([-1 1])
        xlim ([-1 1])
        ylim ([-1.5 1.5]) % set ylim, so that plot frame doesn't change
        yline(out,'--', LineWidth=2)
        text(0.62, -0.25, 'MHWN','FontSize', 14, 'Color', [0.3, 0.3, 0.3])
        xlabel('Domain $x$', Interpreter='latex', FontSize=14)
        ylabel('Elevation $h$', Interpreter='latex', FontSize=14)
        set(gca, "FontSize",14)
    end
    
    % check if front has reached end of domain -> end iteration
    if front_index == 1
        fprintf(['Simulation stopped after %i generations:\nThe ' ...
            'front has reached the end of the domain.\n'], t)
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Main Simulation Figure %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if sum(t == communityPoints) > 0
        figure(10)
        subplotIdx = subplotIdx + 1;
        
        % Plot current densities of all non-extinct species
        subplot(5,4,(13-(subplotIdx*4)):(14-(subplotIdx*4)))
        pl = plot(xx, N, LineWidth=1.5);
        
        % Always use the same colors
        set(pl, {'color'}, num2cell(cmap, 2))
        
        % Set x-axis limits
        xlim([-1 1])
        
        % Add labels and ticks
        if subplotIdx == 1
            xlabelHandle = xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
            xlabelHandle.Position = [0, -0.1588, -1];
        end
        ylabel('$N_{i}$', 'Interpreter', 'latex', 'FontSize', 11)
        
        xticks([-1 -0.5 0 0.5 1])
        xticklabels({'-1', '-0.5', '0', '0.5', '1'})
        yticks([0 0.5 1 1.5])
        yticklabels({'0', '0.5', '1', '1.5'})
        ylim([0 1]) % Set y-axis limits so plot size doesn't change
        
        % Add biomass
        bIdx = bIdx + 1;
        biomass(bIdx, :) = sum(N, 1);
        hold on
        plot(xx, biomass(bIdx, :), 'k', 'LineWidth', 1.5)
        hold off
        
        % Add inset showing timestep
        text(-0.95, 0.95, sprintf('t = %d', t), 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    end
    
    if mod(t,10) == 0 % Every 10 timesteps
        hidx = hidx+1;
        % Save elevation profile
        htime(hidx,:) = h;
    end  
end
disp(t)

%% close video
if visual == true
    fprintf(['Video successfully saved:', videofilename,'.mp4'])
    close(v)
end

%% Travelling waves of elevation profile

figure(10)

% Ensure existing subplot for grayscale is accessed
subplot(5,4,13:20) % Select the subplot where the grayscale plot is needed
ax = gca; % Get the current axes handle
colormap(ax, 'gray'); % Set the colormap for this subplot only

numProfiles = length(hProfilePoints);
cmap = linspace(0.8, 0.2, numProfiles); % Reverse grayscale intensity range (newer = black, older = grey)
pts = round(hProfilePoints / 10);
pts(1) = 1; % First index to 1 / initial elevation profile

% Initialize array for legend labels
legendEntries = cell(length(pts), 1);

cidx = 0;
for hidx = pts
    cidx = cidx + 1;
    
    % Plot each profile
    plot(xx, htime(hidx, :), 'Color', [cmap(cidx) cmap(cidx) cmap(cidx)], 'LineWidth', 2);
    hold on
    
    % Add timestep to legend entries
    legendEntries{cidx} = ['t = ', num2str(hProfilePoints(cidx))];
end
hold off


% Add x and y lines
xline([-1 1])
xlim([-1 1])
ylim([-1 1]) % Set ylim, so that plot frame doesn't change
yline(out, '--', 'LineWidth', 2)

% Add MHWN label
text(0.62, -0.25, 'MHWN', 'FontSize', 14, 'Color', [0.3, 0.3, 0.3])

% Add axis labels
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$h$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, "FontSize", 10)

% Customize ticks
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1', '-0.5', '0', '0.5', '1'})
yticks([-1 0 1])
yticklabels({'-1', '0', '1'})

% Add legend
legend(legendEntries, 'Location', 'northwest', 'FontSize', 8);

%% Add species richness
figure(10)

% Set figure position
set(gcf, 'Position', [90, 30, 620, 550])

subplot(5,4,[3,4,7,8,11,12])
SH = SH(any(~isnan(SH), 2), :); % Remove rows with only NaN values
imagesc(SH);
set(gca, 'YDir', 'normal') % Set axis direction to normal
colormap(jet)

% Adjusting Colorbar
cbar = colorbar;
cbar.Label.String = "D_{l}"; % Add a label to the colorbar
cbar.Label.FontSize = 10; % Smaller font size for the colorbar label
cbar.Label.Rotation = 90; % Rotate the label vertically
cbar.Label.FontName = 'Times New Roman';
cbar.FontSize = 9;

grid off
xlabelHandle = xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 9);
%xlabelHandle.Position = [1000, -60, 0];

% Adjust position of y-axis label
ylabelHandle = ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 9);
%ylabelHandle.Position = [-76, 584, 1]; % Set exact position of the label

xticks([1 n/2 n])
xticklabels({num2str(-L/2), '0', num2str(L/2)})
yticks([1 size(SH,1)])
t_end = t;
yticklabels({num2str(tthres), num2str(t_end)})

set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 12);
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 12);
ax = gca;
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;

% Add position of lower salt marsh boundary to plot:
% Remove NaN values from front
front = front(:,any(~isnan(front), 1));
hold on
% Get gridpoints of each front position
front_gridpoints = (front - min(x)) / dx + 1;
plot(front_gridpoints(tthres+1:end), 0:t_end-tthres, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
hold off

% Adding legend and making it smaller
legend('Marsh edge', 'FontSize', 8); % Smaller font size for the legend

AddLetters2Plots(gcf,{'a)','b)','c)','d)','','e)'},'HShift', -0.05, 'VShift', -0.050)

%%

% safe end time for plots
t_end = t;

% at end of Iteration: transfer all non-extinct species to max_dens matrix
dens_max(N_id,:) = dens_max_iter(:,:);
dens_max_ele(N_id,:) = dens_max_iter_ele(:,:);

% determine elevation at which there is no more sedimentation, but only
% erosion: coresponds to the first value of c, which is <=0
h_thres_index = find(cplot<=0,1,'first');
h_thres = xx(h_thres_index);

% determine number of surviving species
fprintf(['Survived species: ',num2str(size(N,1)), ' out of ', ...
    num2str(m),'\n'])

%% Plot: stress-adjusted growth rates:

figure(4),set(gcf, 'Name', 'stress-adjusted growth rates', ...
    'NumberTitle', 'off');
subplot(1,3,2:3)
plmgr = plot(xx,mgrplot, LineWidth=1.5);
set(plmgr, {'color'}, num2cell(cmapplot,2)) % always use same colours
xlim([-1 1.0])
yline(1,'--',LineWidth=2)
xlabel('$h$', Interpreter='latex', FontSize=12)
ylabel('$r_{i}(x)$', Interpreter='latex', FontSize=12)
ylim([0, Rmaxmax])
%xl4 = xline(xx(h_thres_index),'--', 'HandleVisibility','off', LineWidth=2);
xl5  = xline(0,'--', 'HandleVisibility','off',LineWidth=2);
text(0.05, 3.2, 'MHWN','FontSize', 10, 'Color', [0.3, 0.3, 0.3],'Rotation',90)
text(0.86, 0.3, 'Threshold','FontSize', 10, 'Color', [0.3, 0.3, 0.3],'Rotation',90)
text(1, 0.4, 'elevation','FontSize', 10, 'Color', [0.3, 0.3, 0.3],'Rotation',90)
cbar = colorbar('Direction','reverse');
cbar.Colormap = cmapplot; % use original values of cmapplot
cbar.YTick = [0.05, 0.5, 0.95];
cbar.TickLabels = {'high', 'medium', 'low'};
cbar.Title.String = "Tolerance";
cbar.FontSize=10;
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 10);
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 10);
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})
yticks([0 1 2 3 4])
yticklabels({'0','1','2','3','4'})

% Growth Tolerance trade-off between h0 and Rmax
figure(4),set(gcf, 'Name', 'Trade-off', 'NumberTitle', 'off');
subplot(1,3,1)
plot(hTplot, Rmaxplot, LineWidth=2)
hold on
plot(hT,Rmax,'*', MarkerSize=5, LineWidth=1.5)
hold off
xlabel('$\it{hT}$', 'FontSize',10, Interpreter='latex')
ylabel('Max. growth rate $R$', Interpreter='latex', FontSize=10)
set(gca, "FontSize",10)
xticks([0 0.2 0.4 0.6 0.8])
xticklabels({'0','0.2','0.4','0.6','0.8'})
yticks([1.5 2 2.5 3 3.5 4])
yticklabels({'1.5','2.0','2.5','3.0','3.5','4.0'})
legend('Trade-off curve','Surviving species','Location','southeast')
ylim([1.5 3])

