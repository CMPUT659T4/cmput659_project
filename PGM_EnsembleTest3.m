function Predicted_Y=PGM_EnsembleTest3(Test1,Test2,Test3)
Predicted_Y=zeros(size(Test1));
%disp('gsgafsjafsh');
%size(Predicted_Y,2)
for i=1:length(Predicted_Y)
    %disp('here')
    if Test1(i)==-1 && Test2(i)==-1
        Predicted_Y(i)=-1;
    elseif Test1(i)==1 && Test2(i)==1
        Predicted_Y(i)=1;
    elseif Test1(i)==-1 && Test2(i)==1
        if Test3(i)==-1
            Predicted_Y(i)=1;
        else
            Predicted_Y(i)=-1;
        end
        
    elseif Test1(i)==1 && Test2(i)==-1
        if Test3(i)==-1
            Predicted_Y(i)=1;
        else
            Predicted_Y(i)=-1;
        end
    end
end