function [accu final_res]=combineEnsemble()

[pred1 label]=BestFeaturesPrediction_Mario();
pred2=sgmrf_glasso_mario();
pred3=region_degrees_mario();
pred=pred1+pred2+pred3;
%pred=pred1+pred3;
final_res=zeros(size(pred));
    for i=1:length(pred)
        if pred(i)>0
            final_res(i)=1;
        else
            final_res(i)=-1;
        end
    end
accu=sum(final_res == label) / length(label);
end