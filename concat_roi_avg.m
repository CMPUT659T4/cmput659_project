load('FBIRN/finaldata_AO/SubjectsID_final.mat')
data = zeros(264, 137, 380);
labels = zeros(380, 1);
sites = zeros(380, 1);

roi_dir = dir('ROI_files');
dirs = cell(1, length(roi_dir) - 2);
for folder = 3:length(roi_dir)
    f_name = roi_dir(folder).name;
    dirs{folder - 2} = f_name;
end

for file = 1:length(dirs)
    f_name = dirs{file};
    fprintf([f_name,'\n']);
    
    load(strcat('ROI_files/', f_name));
    data(:, :, file) = ROI_time_series;
    
    s_id = strsplit(f_name, '_');
    s_id = s_id{1};
    site_id = s_id(1:4);
    sites(file) = str2double(site_id);
    
    site_info = getfield(SubjectsID, strcat('SCZ_site', site_id, '_'));
    
    lbl_ind = find(not(cellfun('isempty', strfind(site_info.ID, s_id))));
    labels(file) = site_info.labels(lbl_ind);
end

save(['concat_roi_avg.mat'], 'data', 'labels', 'site')