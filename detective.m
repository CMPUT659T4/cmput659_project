load('concat_roi_avg.mat')

subjects = [];
for s = 1:size(data, 3)
    for r = 1:size(data, 1)
        if squeeze(data(r, :, s)) == zeros(1, size(data, 2))
            subjects = [subjects, r];
            %[s, r]
        end
    end
end
unique([subjects])