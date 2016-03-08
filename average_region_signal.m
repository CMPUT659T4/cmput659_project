function [avg_signal] = average_region_signal(I, indices)
%AVERAGE_REGION_SIGNAL Find the 33 pixal average signal of a region
%   Every region of interest is modeled as sphere of 5 mm radius.
%   This is implemented this as a 3x3 cube + 6 extra voxels
%
%   Authors: Roberto Vega and Neil Borle
    
time_series_counter = 0;
time_series_mat = zeros(33,size(I.img,4));
% First, extract the time series of the 3 x 3 cube
for i = indices(1,1)-1 : indices(1,1) + 1
    for j = indices(2,1)-1 : indices(2,1) + 1
        for k = indices(3,1)-1 : indices(3,1) + 1
            time_series_counter = time_series_counter + 1;
            time_series_mat(time_series_counter,:) = I.img(i,j,k,:);
        end
    end
end

% Extract the time series of the 6 extra voxels
i = indices(1,1);
j = indices(2,1);
k = indices(3,1);
    
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i+2,j,k,:);
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i-2,j,k,:);
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i,j+2,k,:);
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i,j-2,k,:);
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i,j,k+2,:);
time_series_counter = time_series_counter + 1;
time_series_mat(time_series_counter,:) = I.img(i,j,k-2,:);
    
% Get the average of all the signals.
avg_signal = mean(time_series_mat);

end

