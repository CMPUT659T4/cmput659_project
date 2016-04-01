function Accuracy = testCross(data,rho,numFeatures,test_ind)


%numFeatures= 1000;
%y=data(:,end);
%a=find(y==1);
%sch=a(1:38);
%a=find(y==-1);
%hel=a(1:38);
%test_ind=[hel;sch];
train = setdiff([1:size(data,1)],test_ind);
[~,p] = ttest(data(train,:));
[~,I] = sort(p);
ind = I(1:numFeatures);
X = data(train,ind);
%[~,p] = ttest(data(test,:));
%[~,I] = sort(p);
X_ = data(test_ind,ind);
y = data(train,end);
y_ = data(test_ind,end);
method ='varsel_mrf'; %'kataya';%'projected_grad';%'covsel';%'varsel_mrf';
model = MRFC_learn(X, y, method, rho);
[y,pyx] = MRFC_predict(X_, model);

Accuracy = (sum(y==y_)/length(y))*100;
if Accuracy>=70,
    fprintf('Holy Moly Cow!!!!');
    rho
    numFeatures
end;
end