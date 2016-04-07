accu=zeros(20,1);
for i=1:20
accu(i,1)=BestFeaturesPrediction_Final(10);
end
mean(accu)
min(accu)
max(accu)