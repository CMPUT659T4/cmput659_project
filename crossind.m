function ind = crossind(data,numFolds)
L=size(data,1)/4;
indi = crossvalind('Kfold',[1:L],numFolds);
b = repmat(indi,1,4);
b=b';
ind = reshape(b,1,size(b,1)*size(b,2));
ind =ind';
end