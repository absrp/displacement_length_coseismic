%% synthetic fault array to test effect of sampling in 
clear; close all;
% make cracks
% make slip distribution for each crack
% add high-frequency noise (expected for free surface and dynamic mess)
% test effect of segmentation in sampling

rng("default"); % for reproducibility

%% synthetic fracture population
% synthetic profile for one crack test
L_cracks = logspace(log10(1),log10(10000),1000); % meters
scaling = 10^-4; % using scaling factor of 10^-4

% create displacement population for entire crack (tapering from middle,
% triangular) 
interval = 1; % displacements generated every 30 cm
mid_D = L_cracks(500)*scaling; % find maximum displacement for crack of given length 
tip_D = 0; % fix displacement at fault tips at 0
[sampling_points, slip] = triangular_profile(L_cracks(500),tip_D,mid_D,interval);
[rows, columns] = size(slip);
noise = 2 * slip .* (rand(rows, columns)); % add random noise
noisy_slip = slip + noise; % noisy slip profile

figure
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 5, 8]; % [x, y, width, height]

% plot sample displacement profile
subplot(2,1,1)
plot(sampling_points(1:end-1), noisy_slip,'Color',[0.5 0.5 0.5],'linewidth',1.5)
ylabel('Displacement (m)')
xlabel('Distance along the fracture (m)')
set(gca,'FontSize',12)
box on

% create fracture population and estimate max displacement for each crack following scaling
% relationship 
noisy_slip = [];
Lmax_sampled = [];
max_sampled_slip = [];
sampling_points_all = [];
Lsave = [];
thresh = logspace(log10(1),log10(10000),100); % meters

subplot(2,1,2)
for i=1:length(thresh)
    noisy_slip = [];
    Lsave = [];
    threshi = thresh(i);
    for p=1:length(L_cracks)
        if L_cracks(p)>threshi
        interval = 1; % generate displacement estimates every 30 cm 
        mid_D = L_cracks(p)*scaling; % find maximum displacement for crack of given length 
        tip_D = 0;
        [sampling_points, slip] = triangular_profile(L_cracks(i),tip_D,mid_D,interval);
        % add noise
        [rows, columns] = size(slip);
        noise = 2 * slip .* (rand(rows, columns));
        noisy_slipi = slip + noise;
        % save profile
        noisy_slip = [noisy_slip; noisy_slipi];
    else
        continue
        end 
    end 

% generate random indices to simulated segmented mapping
if length(noisy_slip)>1
    numPointsToSample = 10;
    length(noisy_slip)
    randomIndices = randi(length(noisy_slip), 1, numPointsToSample);

% sample points from the array using the random indices
    displ = noisy_slip(randomIndices);

% plot displacement versus length for the randomly sampled segmented cracks
    scatter(repelem(threshi,length(displ)),displ,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
    set(gca,'YScale','log','XScale','log')
    hold on
else
    continue
end
end

% plot relationship for scaling constant 10^-4
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