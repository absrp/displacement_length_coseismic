%%% this script measures the scaling between fracture length and maximum
%%% displacement for four events in the FDHI database
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
     
    % load secondary fracture map for event
    strname = '_secondary_fractures.shp';
    combined_str_sec = append(event,strname);
    lines_secondary = shaperead(combined_str_sec); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% find longest fracture within 10 m of each
%%%%%%%%%%%%%%%%%%%%%%%%%%%% displacement measurement
    % initiate displacement and length arrays
    distance = zeros(length(coords_ref),length(lines_secondary)); 
    sliptest = repmat(slip,1,length(lines_secondary));
    ID = [];
    L = [];

    % measure minimum distance between displacement measurements and every
    % crack, save only instaces where the distance is smaller than 10
    % meters
    for n=1:length(lines_secondary)
        [curvexyx, curvexyy] = wgs2utm(lines_secondary(n).Y,lines_secondary(n).X,11,'N');
        curvexy = [curvexyx' curvexyy'];
        curvexy = rmmissing(curvexy); 
        % nearest neighbor within range
        [~,distance(:,n),~] = distance2curve(curvexy,coords_ref,'linear');
        % turn all numbers where the distance exceeds 10 meters into 0
        distance(distance(:, n) > 10, n) = 0;
        Lline = measure_length(lines_secondary(n).X,lines_secondary(n).Y,11,'N');
        L = [L; Lline];
        ID = [ID; n];
    end

    % repeat L, ID, and coordinate values for the 
    L = repmat(L',length(coords_ref),1);
    ID = repmat(ID',length(coords_ref),1);
    coordsx_percrack = repmat(coordsx,1,length(lines_secondary));
    coordsy_percrack = repmat(coordsy,1,length(lines_secondary));

% % Initialize an array to store the largest values of L for non-zero elements per rows in distance
largest_L = zeros(size(distance, 1), 1);
slip_L =  zeros(size(distance, 1), 1);
ID_L =  zeros(size(distance, 1), 1);
coordsxL = zeros(size(distance, 1), 1);
coordsyL = zeros(size(distance, 1), 1);

% Iterate over each row of distance (distances to all displacement points for
% each crack)
for row = 1:size(distance, 1)
    % find the non-zero elements in the row
    nonzero_indices = find(distance(row, :) ~= 0);

    % for cracks with displacement measurements within 10 meters:
    if ~isempty(nonzero_indices)
        % find the maximum value of L for the non-zero elements in the row
        [largest_L(row), idx] = max(L(row, nonzero_indices));
        idxL = nonzero_indices(idx);
        slip_L(row) = sliptest(row, idxL);
        ID_L(row) = ID(row, idxL);
        coordsxL(row) = coordsx_percrack(row, idxL);
        coordsyL(row) = coordsy_percrack(row, idxL);
    end
end
        select_ID = find(ID_L~=0); % find locations where variables were not zero (i.e. distance < 10 meters)
        largest_L = largest_L(select_ID);
        ID_L = ID_L(select_ID);
        slip_L = slip_L(select_ID);
        coordsxL = coordsxL(select_ID);
        coordsyL = coordsyL(select_ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%% now find locations where multiple slip values have been attributed to the
%%%%%%%%%%%%%%%%%%%%%%%%%%% same fracture and select the maximum slip only 
% find unique values of ind and their corresponding indices
[unique_ind, ~, ind_indices] = unique(ID_L,'stable');

 % find the maximum slip value for each unique index
 max_slip_values = accumarray(ind_indices, slip_L, [], @max);
 max_slip_indices = accumarray(ind_indices, (1:numel(slip_L))', [], @(x) {x(x == max(x))});
 max_slip_indices = cell2mat(max_slip_indices);

 % select related length value for that maximum slip value
 L_maxdisp = largest_L(max_slip_indices); 
 IDselect = ID_L(max_slip_indices);

  % plot D vs L
    figure(3)
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0, 0, 10, 5]; % [x, y, width, height]
    subplot(2,3,i)
    scatter(L_maxdisp, max_slip_values,'MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4)
    hold on
    set(gca,'XScale','log','YScale','log','FontSize',12)
    ylabel('Maximum slip (m)')
    xlabel('Fracture length (m)')

    box on
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
        set(gca,'FontSize',14)
        saveas(gcf,'scaling_maxD_longestL_10m.pdf');


        % test plot to ensure each crack has a single displacemnet
        % measurement and only cracks with displacements are selected
        figure
        
        for b=1:numel(IDselect)
            nind = IDselect(b);
            plot(lines_secondary(nind).X,lines_secondary(nind).Y,'k')
        hold on
        end
        
        scatter(coordsxL(max_slip_indices),coordsyL(max_slip_indices),30,max_slip_values,'filled')
        colormap spring
        colorbar
        axis equal

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
