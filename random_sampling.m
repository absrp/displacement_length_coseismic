%% synthetic fault array to test effect of sampling in 
clear; close all;
%%%%%%%%%%%%%%%%%%%%%%% make cracks
% make slip distribution
% add high-frequency noise (expected for free surface and dynamic mess) 

%% example profile
% create fracture population and estimate max displacement for each crack following scaling
% relationship 
L_cracks = logspace(log10(1),log10(10000),1000); % meters
scaling = 10^-4; % using scaling factor of 10^-4


% create displacement population for entire crack (tapering from middle,
% triangular) Measurements generated every 30 cm. 
interval = 1; % set sampling interval for the crack (default 1/10th of the crack length)
mid_D = L_cracks(500)*scaling; % find maximum displacement for crack of given length 
tip_D = 0;
[sampling_points, slip] = triangular_profile(L_cracks(500),tip_D,mid_D,interval);
[rows, columns] = size(slip);
noise = 2 * slip .* (rand(rows, columns));
noisy_slipi = slip + noise;

figure(3)
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 5, 10]; % [x, y, width, height]

subplot(2,1,1)
plot(sampling_points(1:end-1), noisy_slipi,'Color',[0.5 0.5 0.5],'linewidth',1.5)
ylabel('Displacement (m)')
xlabel('Distance along the fracture (m)')
set(gca,'FontSize',12)
box on

%%

% create fracture population and estimate max displacement for each crack following scaling
% relationship 
L_cracks = logspace(log10(1),log10(10000),1000); % meters
scaling = 10^-4; % using scaling factor of 10^-4
noisy_slip = [];
Lmax_sampled = [];
max_sampled_slip = [];
sampling_points_all = [];
Lsave = [];
thresh = logspace(log10(1),log10(10000),100); % meters
% density 1 m 

subplot(2,1,2)
for i=1:length(thresh)
    noisy_slip = [];
    Lsave = [];
    threshi = thresh(i);
    for p=1:length(L_cracks)
        if L_cracks(p)>threshi
        interval = 0.3; % generate displacement estimates every 30 cm 
        mid_D = L_cracks(p)*scaling; % find maximum displacement for crack of given length 
        tip_D = 0;
        [sampling_points, slip] = triangular_profile(L_cracks(i),tip_D,mid_D,interval);
        % figure
        % plot(sampling_points,slip,'LineWidth',1.5,'Color','k')
        % add noise
        [rows, columns] = size(slip);
        noise = 2 * slip .* (rand(rows, columns));
        noisy_slipi = slip + noise;
        noisy_slip = [noisy_slip; noisy_slipi];
        Lsave = [Lsave; L_cracks(p)];
    else
        continue
        end 
    end 
% Generate random indices
if length(noisy_slip)>1
    numPointsToSample = 10;
    length(noisy_slip)
    randomIndices = randi(length(noisy_slip), 1, numPointsToSample);

% Sample points from the array using the random indices
    displ = noisy_slip(randomIndices);
    scatter(repelem(threshi,length(displ)),displ,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
    set(gca,'YScale','log','XScale','log')
    hold on
else
    continue
end
end

plot(L_cracks,10^-4*L_cracks,'Color',[0.6353    0.0784    0.1843],'LineWidth',1.5)
set(gca,'FontSize',12)
box on
ylabel('Displacement (m)')
xlabel('Fracture length (m)')
saveas(gcf,'scaling_simulation.pdf');

%% function dumpster
function [sampling_points, slip] = triangular_profile(L,tip_D,mid_D,interval)
sampling_points = 0:interval:L; 
% find middle point of sampling points 
% midpt = sampling_points(round(L/2)); % near to middle, by index
% first half
slip_first_half = linspace(tip_D,mid_D,round(length(sampling_points)/2));
slip_second_half = linspace(mid_D,tip_D,round(length(sampling_points)/2)-1);
slip = [slip_first_half'; slip_second_half'];
end 