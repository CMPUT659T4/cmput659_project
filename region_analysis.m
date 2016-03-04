function [] = region_analysis(region_coords)
%REGION_ANALYSIS Generates PGMs based on region coords given
%   This function will go through all the patient files and defined in 
%   Power et al. 2011 and perform clustering. Regions in mmc2.csv
%   This function assumes the following structure of the files:
%   subjectID/session_1/rest_1
%
%   Authors: Neil Borle and Roberto Vega

%dirs = dir();
dirs = {};
row_to_add = [0 0 0 1];

dirs_1 = dir('processed/site_1');
for folder = 3:length(dirs_1)
    f_name = dirs_1(folder).name;
    dirs = [dirs, strcat('processed/site_1/', f_name)];
end

dirs_3 = dir('processed/site_3');
for folder = 3:length(dirs_3)
    f_name = dirs_3(folder).name;
    dirs = [dirs, strcat('processed/site_3/', f_name)];
end

dirs_5 = dir('processed/site_5');
for folder = 3:length(dirs_5)
    f_name = dirs_5(folder).name;
    dirs = [dirs, strcat('processed/site_5/', f_name)];
end

for folder = 1:length(dirs)
    f_name = dirs{folder};
    %f_name = dirs(folder).name;
    %if ~(all(f_name >= '0') && all(f_name <= '9'))
    %    continue;
    %end

    path = strcat(f_name, '/', 'session_1', '/', 'rest_1', '/');
    fprintf([f_name,'\n']);
    I = load_untouch_nii(strcat(path, 'ntswrrest.nii'));
    row_x = I.hdr.hist.srow_x;
    row_y = I.hdr.hist.srow_y;
    row_z = I.hdr.hist.srow_z;
   
    ROI_time_series = zeros(size(region_coords,1), 4);
    % Matrix to convert from indices to MNI coordinates
    conversion_mat = [row_x;row_y;row_z;row_to_add];
    for point = 1:size(region_coords,1)
        indices = floor(conversion_mat\[region_coords(point,:)';1]);
        indices = indices(1:3);
        
        avg_signal = average_region_signal(I, indices);
        cluster_size = cluster(I, indices, avg_signal);
    
        ROI_time_series(point, :) = [cluster_size, indices'];
        fprintf('Finished %i / %i\n', point, size(region_coords,1));
    end

    f_name = strsplit(f_name, '/');
    f_name = f_name(end);
    
    cd ROI_files
    save([f_name{1},'.mat'],'ROI_time_series')
    cd ..
end

