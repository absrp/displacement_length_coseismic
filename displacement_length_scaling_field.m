%% WHAT THIS SCRIPT DOES AND DEPENDENCIES


close all; clear; % clean up before starting

% load displacement data from FDHI database
displacement_data = readtable('data_FDHI.xlsx');
events = {'Landers','EMC', 'HectorMine','Ridgecrest1','Ridgecrest2'}; 

% Find what displacements in the database correlate with each event and 
% save magnitude, x location, and y location per measurement per event 
reflines_all = shaperead('_FDHI_FLATFILE_ECS_rev2.shp'); % ECS lines in FDHI database
%% DISPLACEMENT ANALYSIS

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
     

    % load reference primary fault trace from Rodriguez Padilla and Oskin
    % (2023)
    strname = '_main_rupture.shp';
    combined_str_main = append(event,strname);
    main_rupture = shaperead(combined_str_main); 

    distance = zeros(length(coords_ref),length(main_rupture)); 

for n=1:length(main_rupture)
        [curvexyx, curvexyy] = wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        [xy,distance(:,n),t_a] = distance2curve(curvexy,coords_ref,'linear');
end
 
distance = min(distance,[],2); 
figure(1)
subplot(2,3,i)
% Set the size of the figure
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 10, 6]; % [x, y, width, height]

scatter(distance,slip,'MarkerFaceColor',c,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none') 
ylabel('Slip (m)')
xlabel('Distance from the fault (m)')
set(gca,'YScale','log','XScale','log','FontSize',12)
% Adjust the X and Y axis tick label frequency
xticks('auto'); % Show all X-axis tick labels
yticks('auto'); % Show all Y-axis tick labels

% now attribute each displacement to the longest crack available within 10
% meters, not the nearest
   if i == 4
        title('Ridgecrest foreshock')
    elseif i == 5
        title('Ridgecrest mainshock')
    elseif i == 2
        title('El Mayor-Cucapah')
     elseif i == 3
        title('Hector Mine')     
    else 
        title(event)

    end 
box on
%saveas(gcf,'dispdecay.pdf');


   % locate displacements along the rupture 
   % location of gate along rupture (referenced to ECS files in FDHI database)  
    
    loc_along = [];
    normalized_loc_along = [];
    total_rupturelength = [];


  % find ECS lines for select earthquake
  % find data associated with select earthquake
    EQ_ID = subset_data.EQ_ID(1);
    celllines = struct2cell(reflines_all)'; 
    reflinesloc = find(cell2mat(celllines(:,5)) == EQ_ID); 
    reflines = reflines_all(reflinesloc);
    
    for n = 1:length(coords_ref)
        coordsx = coords_refx(n);
        coordsy = coords_refy(n);
        [total_rupturelengthi,loc_alongi,normalized_loc_alongi] = measure_location_along_rupture_disp(coordsx,coordsy,reflines.X,reflines.Y,11,'N');
        loc_along = [loc_along; loc_alongi];
        normalized_loc_along = [normalized_loc_along; normalized_loc_alongi];
        total_rupturelength = [total_rupturelength; total_rupturelengthi];
    end 

    figure(2)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 5, 10]; % [x, y, width, height]

    subplot(5,1,i);
    scatter(loc_along/1000,slip,20,slip, ...
       'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','None')
    ylabel('Slip (m)')
    xlabel('Distance along the rupture (km)')
    set(gca,'ColorScale','log')
    colormap spring
    colorbar
    xlim([0, max(loc_along)/1000])
    clim([0 5]);
    box on
    saveas(gcf,'slip_profile.pdf');
end 


%% LENGTH ANALYSIS
% magnitude-frequency histograms and spatial

for b=1:length(events)

    event = events{b};
    info = event_info(event);
    c = info{1}; % event color
    
    % subset spreadsheet to event data 
    name = displacement_data.eq_name; % subset event data from FDHI database
    idx = find(strcmp(name,event));
    subset_data = displacement_data(idx,:);

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

    legth_frac_OG = length_frac;

% % break-down fractures into evenly spaced points at 1m spacing and measure
% % distance to nearest point on discretized main rupture
% 
% pt_x = [];
% pt_y = [];
% ID = [];
% 
% distance = [];
% 
% for n=1:numel(lines_secondary)
% %     % generate spline of fracture and resample at 1 m spaced points
%     [pt_x_i,pt_y_i] = subdivide_points(lines_secondary(n).X,lines_secondary(n).Y);
%     pt_x = [pt_x; pt_x_i];
%     pt_y = [pt_y; pt_y_i];
%     ID = [ID; repelem(n,length(pt_x_i))'];
% end
% 
%     for n=1:length(main_rupture)
%         [coords_refx, coords_refy] =  wgs2utm(main_rupture(n).Y,main_rupture(n).X,11,'N');
%         curvexy = [coords_refx' coords_refy'];
%         curvexy = rmmissing(curvexy);
%         [xy,distance(:,n),t_a] = distance2curve(curvexy,[pt_x pt_y],'linear');
%     end
% 
%     dist = min(distance,[],2);
% 
% % find what fractures correspond to what displacement
% 
% if length(dist) ~= length(ID)
%     error('Lengths of distance array and fractures ID array do not match!')
% else
% end
% 
% appended = [dist ID];
% % Extract unique IDs
% uniqueIDs = unique(appended(:, 2));
% % Initialize an array to store the smallest distance per ID
% dist_per_crack = zeros(length(uniqueIDs), 2);
% 
% % Loop through unique IDs and find the smallest distance for each
% for i = 1:length(uniqueIDs)
%     currentID = uniqueIDs(i);
%     % Extract rows with the current ID
%     rows_with_currentID = appended(appended(:, 2) == currentID, :);
%     % Find the smallest distance for the current ID
%     minDist = min(rows_with_currentID(:, 1));
%     % Store the result in the dist_per_crack array
%     dist_per_crack(i, :) = [minDist, currentID];
% end
% 
% figure(3)
% fig = gcf;
% fig.Units = 'inches';
% fig.Position = [0, 0, 10, 5]; % [x, y, width, height]
% subplot(2,3,b)
% scatter(dist_per_crack,length_frac,'MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% hold on
% set(gca,'XScale','log','YScale','log','FontSize',12)
% ylabel('Fracture length (m)')
% xlabel('Distance from fault (m)')
% box on
%    if b == 4
%         title('Ridgecrest foreshock')
%     elseif b == 5
%         title('Ridgecrest mainshock')
%     elseif b == 2
%         title('El Mayor-Cucapah')
%      elseif b == 3
%         title('Hector Mine')     
%     else 
%         title(event)
% 
%    end 
% 
%    saveas(gcf,'length_dist.pdf');
% 
% figure(4)
% fig = gcf;
% fig.Units = 'inches';
% fig.Position = [0, 0, 10, 5]; % [x, y, width, height]
% subplot(2,3,b)
% % bin length data
% % Sample data in length_frac (replace this with your actual data)
% bin_edges = logspace(log10(min(length_frac)), log10(max(length_frac)), 20);
% 
% % Create the histogram counts for the log bins
% hist_counts = histcounts(length_frac, bin_edges);
% 
% histogram('BinCounts',hist_counts,'BinEdges',bin_edges,'FaceColor',c)
% xlabel('Fracture length (m)');
% ylabel('Frequency');
% set(gca,'XScale','log','YScale','log')
% 
% 
% if b == 4
%         title('Ridgecrest foreshock')
%     elseif b == 5
%         title('Ridgecrest mainshock')
%     elseif b == 2
%         title('El Mayor-Cucapah')
%      elseif b == 3
%         title('Hector Mine')     
%     else 
%         title(event)
% 
%    end 
%    saveas(gcf,'length_distribution.pdf');


   % locate fractures along the rupture 
   % location of gate along rupture (referenced to ECS files in FDHI database)  
    
    loc_along = [];
    normalized_loc_along = [];
    total_rupturelength = [];


  % find ECS lines for select earthquake
  % find data associated with select earthquake
    EQ_ID = subset_data.EQ_ID(1);
    celllines = struct2cell(reflines_all)'; 
    reflinesloc = find(cell2mat(celllines(:,5)) == EQ_ID); 
    reflines = reflines_all(reflinesloc);
    
    for n = 1:length(lines_secondary)
        [total_rupturelengthi,loc_alongi,normalized_loc_alongi] = measure_location_along_rupture_frac(lines_secondary(n).X,lines_secondary(n).Y,reflines.X,reflines.Y,11,'N');
        loc_along = [loc_along; loc_alongi];
        normalized_loc_along = [normalized_loc_along; normalized_loc_alongi];
        total_rupturelength = [total_rupturelength; total_rupturelengthi];
        disp(n)
    end 

    figure(5)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 5, 10]; % [x, y, width, height]

    subplot(5,1,b);
    scatter(loc_along/1000,legth_frac_OG,10,legth_frac_OG, ...
       'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','None')
    ylabel('Slip (m)')
    xlabel('Distance along the rupture (km)')
    set(gca,'ColorScale','log')
    colormap spring
    colorbar
    xlim([0, max(loc_along)/1000])
    clim([0 2000]);
    box on
    saveas(gcf,'frac_length_profile.pdf');

end


%% plot maps with secondary fractures, and displacement points
% (color-coded by slip magnitude) 

for i=1:length(events)
    event = events{i};
    c = info{1}; % event color
    
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

% attribute displacement measurements to potential crack hosts

for i=1:length(events)
    event = events{i};
    
    % extract pre-saved color and file name choices for event
    info = event_info(event);
    c = info{1}; % event color

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
    slip_data = [coords_ref slip];

    % load distributed ruptures from FDHI database appendix
    strname = '_secondary_fractures.shp';
    combined_str_sec = append(event,strname);
    lines_secondary = shaperead(combined_str_sec); 

    % get coordinates of fractures (treated as curves here as opposed to
    % discretized points in earlier code block) and find length of each
    % crack. Measure distance between crack and each displacement point
    
    dist = zeros(length(coords_ref),length(lines_secondary)); 
    length_frac = zeros(length(coords_ref),length(lines_secondary)); 
    coords_x = zeros(length(coords_ref),length(lines_secondary)); 
    coords_y = zeros(length(coords_ref),length(lines_secondary)); 

    for n=1:length(lines_secondary)
        [curvexyx, curvexyy] = wgs2utm(lines_secondary(n).Y,lines_secondary(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        [xy,dist(:,n),t_a] = distance2curve(curvexy,coords_ref,'linear'); 
        L = measure_length(lines_secondary(n).X,lines_secondary(n).Y,11,'N');
        length_frac(:,n) = repelem(L,length(xy))';
        coords_x(:,n) = xy(:,1);
        coords_y(:,n) = xy(:,2);
    end

    % first find the minimum distance between a displacement measurement
    % and a crack for all displacement measurements 
    [distance, indices] = min(dist, [], 2);
    coordxtest = [];
    coordytest = [];
    length_frac_test = [];

   for row = 1:size(distance, 1)
       length_frac_testi = length_frac(row, indices(row));
        length_frac_test = [length_frac_test; length_frac_testi'];
    end

    for row = 1:size(distance, 1)
        coordxtesti = coords_x(row, indices(row));
        coordxtest = [coordxtest; coordxtesti'];
    end

    for row = 1:size(distance, 1)
        coordytesti = coords_y(row, indices(row));
        coordytest = [coordytest; coordytesti'];
    end

    % figure
    % scatter(coordxtest,coordytest,30,slip,'filled','MarkerFaceAlpha',0.2)
    % length(coordxtest)
    % length(slip)
    % find overlaping rows between selected coordinates (in order
    % correlating with fractures) and slip data in FDHI database -- still
    % allowing the same crack to be attributed to different displacements
    % (choosing nearest neighbors)
    % xy = [coordxtest coordytest];
    % % % Initialize a third column in xy with NaN
    % xy(:, 3) = NaN;
    % % % Create a logical index for matching rows
    % match_idx = ismember(xy(:, 1:2), slip_data(:, 1:2), 'rows');
    % % % Update the third column of xy with values from slip_data
    % xy(match_idx, 3) = slip_data(match_idx, 3);
    % 
    % figure  
    % hold on
    % scatter(xy(:,1),xy(:,2),30,xy(:,3),'filled','MarkerFaceAlpha',0.2)
    % colormap spring
    % axis equal
    % set(gca,'ColorScale','log','FontSize',14)
    % colorbarLabel = ylabel(colorbar, 'Slip (m)');

    % find what crack length and displacement measurement are associated
    % with those minimum distances
    % find disp related to coords_ref
    
    figure(5)
    subplot(2,3,i)
    scatter(length_frac_test,slip,'MarkerFaceColor',c,'MarkerEdgeColor','None')
    ylabel('Slip (m)')
    xlabel('Fracture length (m)')
    set(gca,'XScale','log','YScale','log')
    
    % calculate scaling relationship D/L w distance
    figure(6)
    subplot(2,3,i)
    scatter(distance,slip/length_frac_test,'MarkerFaceColor',c,'MarkerEdgeColor','None')
    ylabel('D/L')
    xlabel('Distance (m)')
    set(gca,'XScale','log','YScale','log')
    
    % distance array contains the distances between all fractures and all
    % displacement measurements. Let's find displacement measurements
    % within 10 meters of a crack:
    
    
end 


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
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
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

[xy,~,~] = distance2curve(curvexy,[fault_x fault_y],'spline'); % find minimum distance between gate and ECS trace
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
function [total_rupturelength,loc_along,normalized_loc_along] = measure_location_along_rupture_frac(fault_x,fault_y,refline_x,refline_y,zone,hem)

fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));

[coords_gatex, coords_gatey] = subdivide_points(fault_x,fault_y); % discretize fractures into 1 m spaced points for distance measurements
coords_frac = [coords_gatex coords_gatey];

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

spacing = 100; % discretizing rupture into 100 m spaced increments to resample
pt = interparc(0:(spacing/total_rupturelength):1,curvexy_x,curvexy_y,'linear'); 
pt_x = pt(:,1);
pt_y = pt(:,2);
curvexy_dense = [pt_x pt_y];

[xy,distance,~] = distance2curve(curvexy_dense,coords_frac,'linear'); % find minimum distance between gate and ECS trace
[~,index] = min(distance);
% xtest = xy(:,1);
% ytest = xy(:,2); 
% xtest = xtest(index);
% ytest = ytest(index);
% xy = [xtest ytest];

xy = xy(index,:);
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