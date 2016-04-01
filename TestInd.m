function test_ind = TestInd(data,percentage)
    if size(data,1)~=380,
        error('Data format is incompatible with current test format!!');
    end
    numTest = floor((((size(data,1)/4)*percentage)/100)/2);
    test_ind = zeros(numTest*4*2,1);
    testSet = struct('Site',5,'Number_data',5,'Subject',5,'Schizophrenia',5,'Healthy',5,'test_sch',5,'test_hel',5);
    testSet.Site = [1:5];
    testSet.Number_data = [0010,0018,0006,0003,0009];
    testSet.Subject = [23,23,6,21,22];
    testSet.Schizophernia = [10,11,5,10,10];
    testSet.Healthy = [13,12,1,11,12];
    testSet.test_sch = [2,2,1,2,2];
    testSet.test_hel = [3,2,0,2,2];
    yTemp=data(:,end);
    y = zeros(length(yTemp)/4,1);
    for i =1:length(y)
        y(i) = yTemp(i*4);
    end
    pointer = 1;
    count = 0;
    for i=1:5
        ind = find(y(count+1 : count + testSet.Subject(i))== 1);
        ind = ind + count;
        p = randi(length(ind),1,testSet.test_sch(i));
        p = unique(p);
        while (length(p)~=testSet.test_sch(i))
            te = randi(length(ind),1,testSet.test_sch(i)-length(p));
            p = [p te];
            p = unique(p);
        end
        startP =((ind(p)-1).*4)+1;
        startP =[startP;startP+1;startP+2;startP+3];
        ind = sort(startP);
        %ind = reshape(ind,1,size(ind,1)*size(ind,2));
        %ind =ind';
        inii =pointer:pointer+length(ind)-1;
        test_ind(inii) = ind;
        pointer = pointer + length(ind);
        
        ind = find(y(count+1: count + testSet.Subject(i))== -1);
        ind = ind + count;
        count = count + testSet.Subject(i);
        p = randi(length(ind),1,testSet.test_hel(i));
        p = unique(p);
        while (length(p)~=testSet.test_hel(i))
            te = randi(length(ind),1,testSet.test_hel(i)-length(p));
            p = [p te];
            p = unique(p);
        end
        startP =((ind(p)-1).*4)+1;
        startP =[startP;startP+1;startP+2;startP+3];
        ind = sort(startP);
        %ind = reshape(ind,1,size(ind,1)*size(ind,2));
        test_ind(pointer:pointer+length(ind)-1) = ind;
        pointer = pointer + length(ind);
        
    end
    test_ind = sort(test_ind);
end