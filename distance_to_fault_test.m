%%% this script measures the minimum distance between each fracture
%%% measurement and displacement measurement and the principal rupture
%%% trace
% Code requirements:
% Matlab Mapping Toolbox
% Matlab downloadable functions from Mathworks:
% - wsg2utm
% - interparc
% - distance2curve
% Data requirements:
% - Shapefile of secondary fractures
% - Shapefile of primary rupture trace 
% - ECS line from FDHI database (Sarmiento et al., 2021)

close all; clear; % clean up before starting

commandwindow
%% load data 
% load displacement data from FDHI database
displacement_data = readtable('data_FDHI.xlsx');
events = {'Landers','EMC', 'HectorMine','Ridgecrest1','Ridgecrest2'}; 
dist_L = []; % array to store minimum distance between each crack and principal rupture trace
dist_D = []; % array to store minimum distance between each displacement measurement and principal rupture trace

for i=1:length(events)
    event = events{i};
    info = event_info(event);
    c = info{1}; % event color

    % subset spreadsheet to event data 
    name = displacement_data.eq_name; % subset event data from FDHI database
    idx = find(strcmp(name,event));
    subset_data = displacement_data(idx,:);
    type = subset_data.fps_meas_type;
    field = find(strcmp(type,'field')); % select field displacements only for narrow aperture
    subset_data = subset_data(field,:);
    lateral = strcmpi(subset_data.fps_style, {'Right-Lateral'}) | strcmpi(subset_data.fps_style, {'Left-Lateral'});
    subset_data = subset_data(lateral,:);
    slip = subset_data.recommended_net_preferred_for_analysis_meters; % FDHI preferred values
    slipidx = find(slip>0); % avoid artefacts (-999 kinda stuff) 
    slip = slip(slipidx);
    coordsx = subset_data.longitude_degrees(slipidx,:); % get coordinates of measurement
    coordsy = subset_data.latitude_degrees(slipidx,:); 
    EQ_style = subset_data.style;
    measurement_style = subset_data.fps_style;
    [coords_refx, coords_refy] = wgs2utm(coordsy,coordsx,11,'N'); % transform LL to UTM
    coords_ref = [coords_refx coords_refy]; % store coordinates
     
    % load secondary fracture map for event
    strname = '_secondary_fractures.shp';
    combined_str_sec = append(event,strname);
    lines_secondary = shaperead(combined_str_sec); 

    % load reference primary fault trace from Rodriguez Padilla and Oskin
    % (2023)
    strname = '_main_rupture.shp';
    combined_str_main = append(event,strname);
    main_rupture = shaperead(combined_str_main); 


%%%%%%%%%%%%%%%%%%%%%%% measure minimum distance between fractures and
%%%%%%%%%%%%%%%%%%%%%%% rupture

    % initialize variables to populate
    pt_x = [];
    pt_y = [];
    ID = []; % fracture tracker to link segments to full fracture later

    for n=1:numel(lines_secondary)
        % generate spline of fracture and resample at 1 m spaced points
        [pt_x_i,pt_y_i] = subdivide_points(lines_secondary(n).X,lines_secondary(n).Y);
        pt_x = [pt_x; pt_x_i];
        pt_y = [pt_y; pt_y_i];
        ID = [ID; repelem(n,length(pt_x_i))'];
    end

    % initialize distance array
    distance = zeros(length(pt_x),length(main_rupture)); 

    for n=1:length(main_rupture)
        % measure distance between each discretized fracture point and the
        % main rupture
        [coords_refx, coords_refy] =  wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [coords_refx' coords_refy'];
        curvexy = rmmissing(curvexy);
        [xy,distance(:,n),t_a] = distance2curve(curvexy,[pt_x pt_y],'linear');
    end

    dist = min(distance,[],2); % find minimum distance between rupture trace and fracture segment

    % find what fractures correspond to what minimum lengths

    if length(dist) ~= length(ID)
        error('Lengths of distance array and fractures ID array do not match!')
    else
    end
    
    appended = [dist ID];
    % Extract unique IDs
    uniqueIDs = unique(appended(:, 2)); % find fracture IDs for segment groups
    % Initialize an array to store the smallest distance per ID
    dist_per_crack = zeros(length(uniqueIDs),2);
    
    % Loop through unique IDs and find the smallest distance for each
    for i = 1:length(uniqueIDs)
        currentID = uniqueIDs(i);
        % Extract rows with the current ID
        rows_with_currentID = appended(appended(:, 2) == currentID, :);
        % Find the smallest distance for the current ID
        minDist = min(rows_with_currentID(:, 1));
        % Store the result in the dist_per_crack array
        dist_per_crack(i, :) = [minDist,currentID];
    end

    dist_L = [dist_L;dist_per_crack(:,1)];

%%%%%%%%%%%%%%%%%%%%%%% plot distribution of displacements with
%%%%%%%%%%%%%%%%%%%%%%% fault-perpendicular distance

    distance = zeros(length(coords_ref),length(main_rupture)); 

    for n=1:length(main_rupture)
        [curvexyx, curvexyy] = wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        [xy,distance(:,n),t_a] = distance2curve(curvexy,coords_ref,'linear');
    end
 
    distance = min(distance,[],2); 

    dist_D = [dist_D; distance];

end 


%% evaluate distances
close all;

% Log transformation
dist_D_log = log10(dist_D);
dist_L_log = log10(dist_L);

% Compute cumulative sums
[counts_D, edges_D] = histcounts(dist_D_log);
[counts_L, edges_L] = histcounts(dist_L_log);

cum_counts_D = cumsum(counts_D);
cum_counts_L = cumsum(counts_L);

% Calculate the total number of data points
total_points_D = numel(dist_D_log);
total_points_L = numel(dist_L_log);

% Calculate inverse cumulative counts
inv_cum_counts_D = total_points_D - cum_counts_D;
inv_cum_counts_L = total_points_L - cum_counts_L;

figure
subplot(1,2,1)
stairs(edges_D(1:end-1), inv_cum_counts_D, 'Color', [0.5, 0.5, 0.5],'linewidth',2)
ylabel('CDF')
xlabel('Log Distance to displacement measurement (m)')
set(gca, 'YScale', 'log','FontSize',14)

subplot(1,2,2)
stairs(edges_L(1:end-1), inv_cum_counts_L, 'Color', [0.5, 0.5, 0.5],'linewidth',2)
xlabel('Log Distance to fracture measurement (m)')
ylabel('CDF')
set(gca, 'YScale', 'log','FontSize',14)


%% function dumpster
function event_info = event_info(event) 
if strcmp(event,'Landers')
    c = [0.6353    0.0784    0.1843];
    str1 = 'Landers';
    str2 = 'Landers_Scracks.txt';
    str3 = 'Landers_Dpoints.txt';
elseif strcmp(event,'HectorMine')
    c = [0.1647    0.3843    0.2745];
    str1 = 'Hector Mine';
    str2 = 'HectorMine_Scracks.txt';
    str3 = 'HectorMine_Dpoints.txt';
elseif strcmp(event,'Ridgecrest1')
    str1 = 'Ridgecrest foreshock';
    c = [0.8706    0.4902         0];
    str2 = 'Ridgecrest1_Scracks.txt';
    str3 = 'Ridgecrest1_Dpoints.txt';
elseif strcmp(event,'Ridgecrest2')
    str1 = 'Ridgecrest mainshock';
    c = [0.4941    0.1843    0.5569];
    str2 = 'Ridgecrest2_Scracks.txt';
    str3 = 'Ridgecrest2_Dpoints.txt';
elseif strcmp(event,'EMC')
    c = [0    0.6000    0.6000];
    str1 = 'El Mayor Cucapah';
    str2 = 'EMC_Scracks.txt';
    str3 = 'EMC_Dpoints.txt';
else
    slip('Event name does not match database name')
end
    event_info = {c str1 str2 str3};
end 
function [total_rupturelength,loc_along,normalized_loc_along] = measure_location_along_rupture_disp(fault_x,fault_y,refline_x,refline_y,zone,hem)

refline_x = refline_x(~isnan(refline_x));
refline_y = refline_y(~isnan(refline_y));
[curvexy_x, curvexy_y] = wgs2utm(refline_y,refline_x,zone,hem);
curvexy = [curvexy_x' curvexy_y'];

% total length
x_1 = curvexy_x(1:end-1);
x_2 = curvexy_x(2:end);
y_1 = curvexy_y(1:end-1);
y_2 = curvexy_y(2:end);
segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
total_rupturelength = sum(segment);

spacing = 100; % discretizing rupture into 100 m spaced increments
pt = interparc(0:(spacing/total_rupturelength):1,curvexy_x,curvexy_y,'linear'); 
pt_x = pt(:,1);
pt_y = pt(:,2);
curvexy_dense = [pt_x pt_y];

[xy,~,~] = distance2curve(curvexy,[fault_x fault_y],'spline'); % find minimum distance between displacement location and ECS trace
locpt = dsearchn(curvexy_dense,xy);

% segment length
x_1 = pt_x(1:locpt-1);
x_2 = pt_x(2:locpt);
y_1 = pt_y(1:locpt-1);
y_2 =  pt_y(2:locpt);
segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); 
loc_along= sum(segment);

normalized_loc_along = loc_along/total_rupturelength; 
end
function [pt_x,pt_y] = subdivide_points(fault_x,fault_y)
if ~isempty(fault_x) % if statement to deal with empty lines in shapefile
 fault_x = fault_x(~isnan(fault_x));
 fault_y = fault_y(~isnan(fault_y));
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,11,'N');
 %% measure length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
euclidean_matrix = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system
L = sum(euclidean_matrix);
%% create spline of fault 
spacing = 1; % meters
if length(fault_x) > 2
pt = interparc(0:(spacing/L):1,fault_x,fault_y,'linear'); 
pt_x = pt(:,1);
pt_y = pt(:,2);
%spacetest = sqrt((pt_x(1)-pt_x(2)).^2 + (pt_y(1)-pt_y(2)).^2) % check
%that spacing is 1 m 
else
    pt_x = fault_x(1);
    pt_y = fault_y(1);
end
end
end