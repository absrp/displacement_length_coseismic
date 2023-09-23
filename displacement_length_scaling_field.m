close all; clear; % clean up before starting

% load displacement data from FDHI database
displacement_data = readtable('data_FDHI.xlsx');
events = {'Landers','EMC', 'HectorMine','Ridgecrest1','Ridgecrest2'}; 

% Find what displacements in the database correlate with each event and 
% save magnitude, x location, and y location per measurement per event 

%% plot displacement distribution
% spatial
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
    slip = subset_data.recommended_net_preferred_for_analysis_meters; % FDHI preferred values
    slipidx = find(slip>0); % avoid artefacts (-999 kinda stuff) 
    slip = slip(slipidx);
    coordsx = subset_data.longitude_degrees(slipidx,:); % get coordinates of measurement
    coordsy = subset_data.latitude_degrees(slipidx,:); 
    EQ_style = subset_data.style;
    measurement_style = subset_data.fps_style;
    [coords_refx, coords_refy] = wgs2utm(coordsy,coordsx,11,'N'); % transform LL to UTM
    coords_ref = [coords_refx coords_refy]; % store coordinates
     
    % figure(1)
    % subplot(1,5,i)
    % hold on
    % % log bin data 
    % histogram(slip,'FaceColor',c) 
    % ylabel('Frequency')
    % xlabel('Slip (m)')
    % saveas(gcf,'disphist.pdf');
    % 
    if i == 4
        title('Ridgecrest mainshock')
    elseif i == 5
        title('Ridgecrest foreshock')
    else 
        title(event)

    end 

    if i == 5
        p == 4;
    else 
        p = i;
    end

    % load reference primary fault trace from Rodriguez Padilla and Oskin
    % (2023)
    strname = '_main_rupture.shp';
    combined_str_main = append(event,strname);
    main_rupture = shaperead(combined_str_main); 

    figure(1)
    subplot(2,2,p)
    scatter(coordsx,coordsy,20,'MarkerFaceColor',c,'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none') 
    hold on
    for n=1:length(main_rupture)
        plot(main_rupture(n).X,main_rupture(n).Y,'k','linewidth',1.5)
        hold on
        %set(gca,'FontSize',14)
    end
    axis equal
    % box on
    xlabel('Lon')
    ylabel('Lat')

    if p == 1
        text(-116.95, 34.68, 'Landers', 'VerticalAlignment', 'top')
    elseif p == 2
       text(-115.5, 32.715, 'El Mayor-Cucapah', 'VerticalAlignment', 'top')
    elseif p == 3
       text(-116.2, 34.76, 'Hector Mine', 'VerticalAlignment', 'top')
    elseif p == 4
       text(-117.59, 35.8, 'Mainshock', 'VerticalAlignment', 'top','Color', [0.4941    0.1843    0.5569])
       text(-117.7, 35.659, 'Foreshock', 'VerticalAlignment', 'top','Color', [0.8706    0.4902         0])
       text(-117.46, 35.92, 'Ridgecrest', 'VerticalAlignment', 'top')
    end

%saveas(gcf,'dispmap.pdf');

distance = zeros(length(coords_ref),length(main_rupture)); 


for n=1:length(main_rupture)
        [curvexyx, curvexyy] = wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        [xy,distance(:,n),t_a] = distance2curve(curvexy,coords_ref,'linear');
end
 
distance = min(distance,[],2); 
figure(2)
subplot(2,3,i)
scatter(distance,slip,'MarkerFaceColor',c,'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none') 
ylabel('Slip (m)')
xlabel('Distance from fault (m)')
set(gca,'YScale','log','XScale','log')

%saveas(gcf,'dispdecay.pdf');
end 

%% plot length distribution
% histograms and spatial

for b=1:length(events)

    event = events{b};
    info = event_info(event);
    c = info{1}; % event color

    % load reference primary fault trace from Rodriguez Padilla and Oskin
    % (2023)
    strname = '_main_rupture.shp';
    combined_str_main = append(event,strname);
    main_rupture = shaperead(combined_str_main); 

    % load distributed ruptures from FDHI database appendix
    strname = '_secondary_fractures.shp';
    combined_str_sec = append(event,strname);
    lines_secondary = shaperead(combined_str_sec); 

    length_frac = [];
    % measure lenght of fractures

    for n=1:length(lines_secondary)
        [L] = measure_length(lines_secondary(n).X,lines_secondary(n).Y,11,'N');
        length_frac = [length_frac; L'];
    end 


% break-down fractures into evenly spaced points at 1m spacing and measure
% distance to nearest point on discretized main rupture

pt_x = [];
pt_y = [];
ID = [];

distance = [];

for n=1:numel(lines_secondary)
%     % generate spline of fracture and resample at 1 m spaced points
    [pt_x_i,pt_y_i] = subdivide_points(lines_secondary(n).X,lines_secondary(n).Y);
    pt_x = [pt_x; pt_x_i];
    pt_y = [pt_y; pt_y_i];
    ID = [ID; repelem(n,length(pt_x_i))'];
end

    for n=1:length(main_rupture)
        [coords_refx, coords_refy] =  wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [coords_refx' coords_refy'];
        curvexy = rmmissing(curvexy);
        [xy,distance(:,n),t_a] = distance2curve(curvexy,[pt_x pt_y],'linear');
    end
    
    dist = min(distance,[],2);

%% find what fractures correspond to what displacement

if length(dist) ~= length(ID)
    error('Lengths of distance array and fractures ID array do not match!')
else
end

appended = [dist ID];
% Extract unique IDs
uniqueIDs = unique(appended(:, 2));
% Initialize an array to store the smallest distance per ID
dist_per_crack = zeros(length(uniqueIDs), 2);

% Loop through unique IDs and find the smallest distance for each
for i = 1:length(uniqueIDs)
    currentID = uniqueIDs(i);
    % Extract rows with the current ID
    rows_with_currentID = appended(appended(:, 2) == currentID, :);
    % Find the smallest distance for the current ID
    minDist = min(rows_with_currentID(:, 1));
    % Store the result in the dist_per_crack array
    dist_per_crack(i, :) = [minDist, currentID];
end

figure(3)
subplot(2,3,b)
scatter(dist_per_crack,length_frac,'MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
hold on
set(gca,'XScale','log','YScale','log')
ylabel('Fracture length (m)')
xlabel('Distance from fault (m)')

%saveas(gcf,'length_dist.pdf');

end


%% plot maps with secondary fractures, and displacement points
% (color-coded by slip magnitude) 

for i=1:length(events)
    event = events{i};
    
    % extract pre-saved color and file name choices for event
    info = event_info(event);

    % subset spreadsheet to event data 
    name = displacement_data.eq_name; 
    idx = find(strcmp(name,event));
    subset_data = displacement_data(idx,:);
    type = subset_data.fps_meas_type;
    field = find(strcmp(type,'field'));
    subset_data = subset_data(field,:);
    slip = subset_data.recommended_net_preferred_for_analysis_meters;
    slipidx = find(slip>0); 
    slip = slip(slipidx);
    coordsx = subset_data.longitude_degrees(slipidx,:);
    coordsy = subset_data.latitude_degrees(slipidx,:); 
    EQ_style = subset_data.style;
    measurement_style = subset_data.fps_style;
    [coords_refx, coords_refy] = wgs2utm(coordsy,coordsx,11,'N');
    coords_ref = [coords_refx coords_refy];
  

    % load distributed ruptures from FDHI database appendix
    strname = '_secondary_fractures.shp';
    combined_str_sec = append(event,strname);
    lines_secondary = shaperead(combined_str_sec); 

    figure
    hold on 

    % plot fractures for event 
    for n=1:length(lines_secondary)
        plot(lines_secondary(n).X,lines_secondary(n).Y,'k')
    end

    scatter(coordsx,coordsy,30,slip,'filled','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none')
    colorbar
    set(gca,'ColorScale','log','FontSize',14)
    colorbarLabel = ylabel(colorbar, 'Slip (m)');

    xlabel('Lon')
    ylabel('Lat')
    colormap spring
    axis equal

     if i == 5
        title('Ridgecrest mainshock')
     elseif i == 2
         title('El Mayor-Cucapah')
    elseif i == 4
        title('Ridgecrest foreshock')
    else 
        title(event)

     end 

end 


%% scaling


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

function [fracture_length] = measure_length(fault_x,fault_y,zone,hem)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,11,'N');
% calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
fracture_length = sum(segment_length);
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

