mode = 'regions'; %'regions';

if strcmp(mode, 'regions')
    M = csvread('mmc2.csv');
    region_analysis(M(:, 2:end))
elseif strcmp(mode, 'degrees')
    % insert degrees code
end