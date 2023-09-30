%% synthetic fault array to test effect of sampling in 
clear; close all;
%%%%%%%%%%%%%%%%%%%%%%% make cracks
% make slip distribution
% add high-frequency noise (expected for free surface and dynamic mess) 
% allow scaling to fluctuate?

% create fracture population and estimate max displacement for each crack following scaling
% relationship 
L_cracks = logspace(log10(1),log10(10000),100); % meters
scaling = 10^-4; % using scaling factor of 10^-4
noisy_slip = [];
Lmax_sampled = [];
max_sampled_slip = [];
sampling_points_all = [];

% create displacement population for entire crack (tapering from middle,
% triangular) Measurements generated every 10 cm. 

for i=1:length(L_cracks)
    interval = L_cracks(i)*0.01; % set sampling interval for the crack (default 1/10th of the crack length)
    mid_D = L_cracks(i)*scaling; % find maximum displacement for crack of given length 
    tip_D = 0;
    [sampling_points, slip] = triangular_profile(L_cracks(i),tip_D,mid_D,interval);
    % figure
    % plot(sampling_points,slip,'LineWidth',1.5,'Color','k')
    % add noise
    [rows, columns] = size(slip);
    noise = 2 * slip .* (rand(rows, columns));
    noisy_slipi = slip + noise;
    noisy_slip = [noisy_slip; noisy_slipi'];
    sampling_points_all = [sampling_points_all; sampling_points];
end

figure(3)
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 5, 10]; % [x, y, width, height]

subplot(3,1,1)
plot(sampling_points(end,:), noisy_slip(end,:),'Color',[0.5 0.5 0.5],'linewidth',1.5)
ylabel('Slip (m)')
xlabel('Distance along the fracture (m)')
title('Sample synthetic profile with random noise')
set(gca,'FontSize',12)
box on

%%%%%%%%%%%%%%%%%%%%%%%%% random sampling
max_points_to_sample = 5;  % upper limit of possible sampling


for i = 1:length(L_cracks)
    % Sample a row from noisy_slip
    sampled_row = noisy_slip(i, :);
    
    % Randomly sample a number of points (at least one, up to max_points_to_sample)
    num_points_to_sample = randi([1, max_points_to_sample]);
    sampled_points = datasample(sampled_row, num_points_to_sample, 'Replace', false);
    
    % Calculate the maximum value in the sampled_points array
    max_sampled_point = max(sampled_points);   
    max_sampled_slip = [max_sampled_slip; max_sampled_point];

end

subplot(3,1,2)
scatter(L_cracks,max_sampled_slip,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none')
hold on
plot(L_cracks,10^-4*L_cracks,'Color',[0.6353    0.0784    0.1843],'LineWidth',1.5)
ylabel('D_{max}(m)')
xlabel('Length (m)')
set(gca,'YScale','log','XScale','log','FontSize',12)
title('Random sampling of displacement along fracture')
xlim([5*10^-3 10^4])
ylim([10^-6 2.5])
box on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break down fracture into sub-fractures
max_sampled_slip = [];
sampled_length_segment = [];
num_points_to_sample = 2;
max_points_to_sample = 1;

for i = 1:length(L_cracks)
    % Sample the row i from the sampling_points_all array
    sampled_row = sampling_points_all(i, :);
    tosample = 1:length(sampled_row);
    % Randomly sample two points from the sampled row
    sampled_points = datasample(tosample, num_points_to_sample, 'Replace', false);
    
    % Sort the sampled points to ensure the correct order
    sampled_points = sort(sampled_points);
    
    % Crop the sampling_points_all array between the two sampled points
    fracture_segment = sampled_row(sampled_points(1):sampled_points(2));

    % measure length of fracture segment 
    sampled_length_segmenti = fracture_segment(2)-fracture_segment(1);
    sampled_length_segment = [sampled_length_segment; sampled_length_segmenti];
   
    % Crop the noisy_slip array between the same two sampled points
    noisy_slip_segment = noisy_slip(i, sampled_points(1):sampled_points(2));
    
    % now sample a few points and keep the maximum slip for that segment
    num_points_to_sample_slip = randi([1, max_points_to_sample]);
    slip_sampled_points = datasample(noisy_slip_segment, num_points_to_sample_slip, 'Replace', false);
    
    % Calculate the maximum value in the sampled_points array
    max_noisy_slip_segment = max(slip_sampled_points);
    max_sampled_slip = [max_sampled_slip; max_noisy_slip_segment];

end

subplot(3,1,3)
scatter(sampled_length_segment,max_sampled_slip,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none')
hold on
xsegments = linspace(10^-1,max(sampled_length_segment),10000);
plot(xsegments,10^-4*xsegments,'Color',[0.6353    0.0784    0.1843],'LineWidth',1.5)
ylabel('D_{max}(m)')
xlabel('Length (m)')
xlim([5*10^-3 10^4])
ylim([10^-6 2.5])
set(gca,'YScale','log','XScale','log','FontSize',12)
box on
title('Breaking down of fracture into segments')
saveas(gcf,'scaling_simulation.pdf');

% add minimum displacement measured 1 mm 
% add shortest crack measured 


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