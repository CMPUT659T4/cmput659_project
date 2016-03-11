clear
load('ebc_subj14_preproc.mat'); 
    NoOfVoxels = size(fmridata,1);
    Length = size(fmridata,2);
 
k =22; %22;
thresh  = 0.8; %feature-fmri correlation
%c=zeros(25)
%for k=1:25
 % Finding correction of feature vector for each voxel
        disp('Finding correlation');
        for ct=1:NoOfVoxels
            xy = corrcoef(featureVector(k,1:704),fmridata(ct,1:704));
            corrs(ct) = xy(1, 2);
        end
%c(k)=max(abs(corrs))
%end
    % Choosing voxels which is highly correlated to feature vector
        disp('Choosing voxels which is highly correlated to feature vector');
    	[bog,inds] = find(abs(corrs) > thresh*max(abs(corrs)));
        disp(sprintf('No of voxels choosen for correlation is %i %i',size(inds,2)));


%%%%

time_points=704;
train_data13=zeros(1+size(inds,2),time_points);
train_data13(1,1:time_points)= featureVector(k,1:time_points);
train_data13(2:(size(inds,2)+1),1:time_points) = fmridata(inds,1:time_points);

m=transpose(mean(train_data13'));
for j=1:(size(inds,2)+1),
    for i=1:time_points,
        train_data13(j,i)= train_data13(j,i) - m(j) ;
    end
end

%m=transpose(var(train_data13'));
%for j=1:(size(inds,2)+1),
%    for i=1:time_points,
%        train_data13(j,i)= train_data13(j,i)/sqrt(m(j)) ;
%    end
%end




data= train_data13'; %% to scale down
save('c2007_22_th06.txt','data','-ASCII','-TABS');


test_points=1408;
test_data13=zeros(1+size(inds,2),time_points);
test_data13(1,1:time_points)= featureVector(k,(time_points+1):test_points);
test_data13(2:(size(inds,2)+1),1:time_points) = fmridata(inds,(time_points+1):test_points);

m=transpose(mean(test_data13'));
for j=1:(size(inds,2)+1),
    for i=1:time_points,
        test_data13(j,i)= test_data13(j,i) - m(j) ;
    end
end

%m=transpose(var(test_data13'));
%for j=1:(size(inds,2)+1),
%    for i=1:time_points,
%        test_data13(j,i)= test_data13(j,i) /sqrt(m(j)) ;
%    end
%end

test_data=test_data13'; %% to scale down
save('c2007_22_th06_test.txt','test_data','-ASCII','-TABS');

 disp(sprintf('finished'));




%%%% Smoothing:


