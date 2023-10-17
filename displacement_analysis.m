%%% this script analyzes the distribution of displacements measured in the
%%% field and stored in the FDHI database
% Code requirements:
% Matlab Mapping Toolbox
% Matlab downloadable functions from Mathworks:
% - wsg2utm
% - interparc
% - distance2curve
% Data requirements:
% - Shapefile of primary rupture trace 
% - ECS line from FDHI database (Sarmiento et al., 2021)

close all; clear; % clean up before starting

%% load data 
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


%%%%%%%%%%%%%%%%%%%%%%% plot distribution of displacement magnitudes along
%%%%%%%%%%%%%%%%%%%%%%% strike of the rupture
    
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
    %saveas(gcf,'slip_profile.pdf');
    
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
