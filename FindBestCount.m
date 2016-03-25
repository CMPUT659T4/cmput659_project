
lim=100;
accuracy=zeros(1,lim);
wrongHeal=zeros(1,lim);
wrongPati=zeros(1,lim);
for feat=1:lim
    fprintf('Calculating the value for %d variables',feat);
    [a b c Met WrongMAt GoodMat]=predictBasedOnBestFeaturesInsight(feat,200);
    accuracy(feat)=mean(a);
    wrongHeal(feat)=mean(b);
    wrongPati(feat)=mean(c);
end