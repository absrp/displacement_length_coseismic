%%% this script analyzes the length-magnitude distribution of fractures
%%% from surface rupture maps and the spatial distribution of fracture
%%% length. 
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

%% load data 
% load displacement data from FDHI database
displacement_data = readtable('data_FDHI.xlsx');
events = {'Landers','EMC', 'HectorMine','Ridgecrest1','Ridgecrest2'}; 

% Find what displacements in the database correlate with each event and 
% save magnitude, x location, and y location per measurement per event 
reflines_all = shaperead('_FDHI_FLATFILE_ECS_rev2.shp'); % ECS lines in FDHI database

%% length distribution analysis 
for b=1:length(events)

    % extract event information from FDHI info database 
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
% 
% %%%%%%%%%%%%%%%%%%%%%%% measure fracture length

    length_frac = [];

    for n=1:length(lines_secondary)
        [L] = measure_length(lines_secondary(n).X,lines_secondary(n).Y,11,'N');
        length_frac = [length_frac; L'];
    end 


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

    figure(3)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 10, 5]; % [x, y, width, height]
    subplot(2,3,b)
    scatter(dist_per_crack(:,1),length_frac,'MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
    hold on
    set(gca,'XScale','log','YScale','log','FontSize',14)
    ylabel('Fracture length (m)')
    xlabel('Distance from fault (m)')

box on
    
       %saveas(gcf,'length_dist.pdf');

%%%%%%%%%%%%%%%%%%%%%%% plot length-magnitude distribution of fractures
        box on
       if b == 4
            title('Ridgecrest foreshock')
            Lmin = 1.8;
        elseif b == 5
            title('Ridgecrest mainshock')
            Lmin = 2;
        elseif b == 2
            title('El Mayor-Cucapah')
            Lmin = 20;
         elseif b == 3
            title('Hector Mine')   
            Lmin = 10;
        else 
            title(event)
            Lmin = 60;
       end 
    figure(2)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 10, 5]; % [x, y, width, height]
    subplot(2,3,b)
    % bin length data
    bin_edges = logspace(log10(min(length_frac)), log10(max(length_frac)), 20);
    % Create the histogram counts for the log bins
    hist_counts = histcounts(length_frac, bin_edges);
    % plot hisgram
    histogram('BinCounts',hist_counts,'BinEdges',bin_edges,'FaceColor',c,'FaceAlpha',0.6)
    hold on
    xline(Lmin,'Color',c,'LineWidth',2)
    xlabel('Fracture length (m)');
    ylabel('Frequency');
    set(gca,'XScale','log','YScale','log','FontSize',14)

    if b == 4
            title('Ridgecrest foreshock')
        elseif b == 5
            title('Ridgecrest mainshock')
        elseif b == 2
            title('El Mayor-Cucapah')
         elseif b == 3
            title('Hector Mine')     
        else 
            title(event)
       end 
       saveas(gcf,'length_distribution.pdf');

%%%%%%%%%%%%%%%%%%%%%%% plot distribution of fracture lengths along
%%%%%%%%%%%%%%%%%%%%%%% strike of the rupture

    % initialize variables
    loc_along = [];
    normalized_loc_along = [];
    total_rupturelength = [];

    % find ECS lines for select earthquake
    % find data associated with select earthquake
    EQ_ID = subset_data.EQ_ID(1);
    celllines = struct2cell(reflines_all)'; 
    reflinesloc = find(cell2mat(celllines(:,5)) == EQ_ID); 
    reflines = reflines_all(reflinesloc);
    refline_x = reflines.X; 
    refline_y = reflines.Y; 
    refline_x = refline_x(~isnan(refline_x));
    refline_y = refline_y(~isnan(refline_y));   
    [curvexy_x, curvexy_y] = wgs2utm(refline_y,refline_x,11,'N');
    curvexy = [curvexy_x' curvexy_y'];

    % measure total length
    x_1 = curvexy_x(1:end-1);
    x_2 = curvexy_x(2:end);
    y_1 = curvexy_y(1:end-1);
    y_2 = curvexy_y(2:end);
    segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
    total_rupturelength = sum(segment);

    % densify ECS line
    spacing = 1; % discretizing rupture into 100 m spaced increments to resample
    pt = interparc(0:(spacing/total_rupturelength):1,curvexy_x,curvexy_y,'linear'); 
    pt_x = pt(:,1);
    pt_y = pt(:,2);

    % check whether ECS goes NS, SN, EW, or WE for each event
    if pt_y(1)>pt_y(end)
        disp(event)
        disp('NS')
    else
        disp(event)
        disp('SN')
    end

    if pt_x(end)>pt_x(1)
        disp(event)
        disp('WE')
    else
        disp(event)
        disp('EW')
    end


    curvexy_dense = [pt_x pt_y];
    length_frac = [];

    % find closest location between fracture and ECS line for each fracture
    for n = 1:length(lines_secondary)
        [loc_alongi,normalized_loc_alongi] = measure_location_along_rupture_frac(lines_secondary(n).X,lines_secondary(n).Y,curvexy_dense,total_rupturelength);
        loc_along = [loc_along; loc_alongi];
        normalized_loc_along = [normalized_loc_along; normalized_loc_alongi];
        [L] = measure_length(lines_secondary(n).X,lines_secondary(n).Y,11,'N');
        length_frac = [length_frac; L'];
    end 

    % plot distribution of fractures along strike of the rupture
    figure(3)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 5, 13]; % [x, y, width, height]

    subplot(5,1,b);
    scatter(loc_along/1000,length_frac,10,length_frac, ...
       'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','None')
    ylabel('Length (m)')
    xlabel('Distance along the rupture (km)')
    set(gca,'ColorScale','log')
    set(gca,'YScale','log')
    colormap spring
    colorbar
    xlim([0, max(loc_along)/1000])
    ylim([0.1, 3*10^3])
    size(get(gca, 'YTick'))
    set(gca,'ColorScale','log','FontSize',16)
    clim([0.1 1500]);
    box on
    saveas(gcf,'frac_length_profile.pdf');

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
function [loc_along,normalized_loc_along] = measure_location_along_rupture_frac(fault_x,fault_y,curvexy_dense,total_rupturelength)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));

[coords_gatex, coords_gatey] = subdivide_points(fault_x,fault_y); % discretize fractures into 1 m spaced points for distance measurements
coords_frac = [coords_gatex coords_gatey];

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
pt_x = curvexy_dense(:,1);
pt_y = curvexy_dense(:,2);
x_1 = pt_x(1:locpt-1);
x_2 = pt_x(2:locpt);
y_1 = pt_y(1:locpt-1);
y_2 =  pt_y(2:locpt);
segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); 
loc_along= sum(segment);

normalized_loc_along = loc_along/total_rupturelength; 
end