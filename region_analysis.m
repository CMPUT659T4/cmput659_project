function [] = region_analysis(region_coords)
%REGION_ANALYSIS Generates PGMs based on region coords given
%   This function will go through all the patient files and defined in 
%   Power et al. 2011 and perform clustering. Regions in mmc2.csv
%   This function assumes the following structure of the files:
%   subjectID/session_1/rest_1
%
%   Authors: Neil Borle and Roberto Vega

%dirs = dir();
row_to_add = [0 0 0 1];

lowres_dir = dir('FBIRN/finaldata_AO/lowresRef');
dirs = cell(1, length(lowres_dir) - 3);
for folder = 3:length(lowres_dir)-1
    f_name = lowres_dir(folder).name;
    dirs{folder - 2} = strcat('FBIRN/finaldata_AO/lowresRef/', f_name);
end

for folder = 1:length(dirs)
    f_name = dirs{folder};
    fprintf([f_name,'\n']);
    
    I = load_untouch_nii(f_name);
    row_x = I.hdr.hist.srow_x;
    row_y = I.hdr.hist.srow_y;
    row_z = I.hdr.hist.srow_z;
    dim = size(I.img);
   
    ROI_time_series = zeros(size(region_coords,1), dim(end));
    % Matrix to convert from indices to MNI coordinates
    conversion_mat = [row_x;row_y;row_z;row_to_add];
    for point = 1:size(region_coords,1)
        indices = floor(conversion_mat\[region_coords(point,:)';1]);
        indices = indices(1:3);
        
        avg_signal = average_region_signal(I, indices);
        ROI_time_series(point, :) = avg_signal;
        fprintf('Finished %i / %i\n', point, size(region_coords,1));
    end

    f_name = strsplit(f_name, '/');
    f_name = strsplit(f_name{end}, '.');
    f_name = f_name{1};
    
    cd ROI_files
    save([f_name,'.mat'],'ROI_time_series')
    cd ..
end

