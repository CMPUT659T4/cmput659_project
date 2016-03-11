%%%%%%%%%%%%%% scale-free network results
%%% plot lambda vs b

data1=load('SF/res_sum_table_exp_scalar_SF.txt');
i30 = find(data1(:,1)==30);
i50 = find(data1(:,1)==50);
i500 = find(data1(:,1)==500);
i1000 = find(data1(:,1)==1000);
figure(1); 
loglog(data1(i30(1:5),3),data1(i30(1:5),5),'bx-');
hold on; 
loglog(data1(i50(1:5),3),data1(i30(1:5),5),'gv-');
loglog(data1(i500(1:5),3),data1(i500(1:5),5),'k.-');
loglog(data1(i1000(1:5),3),data1(i1000(1:5),5),'ro-');

xlabel('b','FontSize', 16); ylabel('lambda','FontSize', 16);
legend('N=30','N=50','N=500','N=1000');
title('scale-free networks','FontSize', 16);



%%%% plot ROC curves
data=load('SF/res_sum_table_ROC_SF.txt');
dataCV=load('SF/res_sum_table_CV_SF.txt');
dataB=load('SF/res_sum_table_banerjee_SF.txt');
dataV=load('SF/res_sum_table_exp_vector_SF.txt');

i30 = find(dataCV(:,1)==30);
i50 = find(dataCV(:,1)==50);
i500 = find(dataCV(:,1)==500);
i1000 = find(dataCV(:,1)==1000);
 
 
i30sp1 = find(data(:,1)==30 & data(:,2) == 0.0618);
i30sp2 = find(data(:,1)==30 & data(:,2) == 0.2194);
i30sp3 = find(data(:,1)==30 & data(:,2) == 0.3138);
i50sp1 = find(data(:,1)==50 & data(:,2) == 0.0618);
i50sp2 = find(data(:,1)==50 & data(:,2) == 0.2194);
i50sp3 = find(data(:,1)==50 & data(:,2) == 0.3138);
i500sp1 = find(data(:,1)==500 & data(:,2) == 0.0618);
i500sp2 = find(data(:,1)==500 & data(:,2) == 0.2194);
i500sp3 = find(data(:,1)==500 & data(:,2) == 0.3138);
i1000sp1 = find(data(:,1)==1000 & data(:,2) == 0.0618);
i1000sp2 = find(data(:,1)==1000 & data(:,2) == 0.2194);
i1000sp3 = find(data(:,1)==1000 & data(:,2) == 0.3138);

%30 0.0618 0.1 1.510227 0.6287225 0.6066742
%30 0.2194 0.1 1.537929 0.6191124 0.6034565
%30 0.3138 0.1 1.541404 0.6189902 0.6055202

figure(2);
hold on;
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30sp1,6),data(i30sp1,5),'k-.');
plot(dataCV(i30(1),6),dataCV(i30(1),5),'kx','MarkerSize',12);
plot(dataB(i30(1),6),dataB(i30(1),5),'kd','MarkerSize',12);
%plot(dataV(i30(1),6),dataV(i30(1),5),'kv','MarkerSize',12);


plot(data(i50sp1,6),data(i50sp1,5),'g-.');
plot(dataCV(i50(1),6),dataCV(i50(1),5),'gx','MarkerSize',12);
plot(dataB(i50(1),6),dataB(i50(1),5),'gd','MarkerSize',12);
%plot(dataV(i50(1),6),dataV(i50(1),5),'gv','MarkerSize',12);
 

plot(data(i500sp1,6),data(i500sp1,5),'b--');
plot(dataCV(i500(1),6),dataCV(i500(1),5),'bx','MarkerSize',12);
plot(dataB(i500(1),6),dataB(i500(1),5),'bd','MarkerSize',12);
%plot(dataV(i500(1),6),dataV(i500(1),5),'bv','MarkerSize',12);

plot(data(i1000sp1,6),data(i1000sp1,5),'r-');
plot(dataCV(i1000(1),6),dataCV(i1000(1),5),'rx','MarkerSize',12);
plot(dataB(i1000(1),6),dataB(i1000(1),5),'rd','MarkerSize',12);
%plot(dataV(i1000(1),6),dataV(i1000(1),5),'rv','MarkerSize',12);
title('scale-free networks','FontSize', 16);
%legend('N=30','CV(N=30)','Banerjee(N=30)','Vector(N=30)','N=50','CV(N=50)','Banerjee(N=50)','Vector(N=50)','N=500','CV(N=500)','Banerjee(N=500)','Vector(N=500)','N=1000','CV(N=1000)','Banerjee(N=1000)','Vector(N=1000)');
legend('N=30','CV(N=30)','Banerjee(N=30)','N=50','CV(N=50)','Banerjee(N=50)','N=500','CV(N=500)','Banerjee(N=500)','N=1000','CV(N=1000)','Banerjee(N=1000)');


figure(3);
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30sp2,6),data(i30sp2,5),'k-.');
hold on;
plot(data(i500sp2,6),data(i500sp2,5),'b--');
plot(data(i1000sp2,6),data(i1000sp2,5),'r-');

figure(4);
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30sp3,6),data(i30sp3,5),'k-.');
hold on;
plot(data(i500sp3,6),data(i500sp3,5),'b--');
plot(data(i1000sp3,6),data(i1000sp3,5),'r-');

%%%%%%%%%%%%%% random network results

data=load('random/res_sum_table_exp_scalar_random.txt');
i30 = find(data(:,1)==30);
i50 = find(data(:,1)==50);
i500 = find(data(:,1)==500);
i1000 = find(data(:,1)==1000);
 
figure(11); 
loglog(data(i30(1:5),3),data(i30(1:5),5),'bx-');
hold on; 
loglog(data(i50(1:5),3),data(i30(1:5),5),'gv-');
loglog(data(i500(1:5),3),data(i500(1:5),5),'k.-');
loglog(data(i1000(1:5),3),data(i1000(1:5),5),'ro-');

xlabel('b','FontSize', 16); ylabel('lambda','FontSize', 16);
legend('N=30','N=50','N=500','N=1000');
title('random networks','FontSize', 16);

%%%% plot ROC curves
%data=load('SF/res_sum_table_ROC_SF.txt');
%dataCV=load('SF/res_sum_table_CV_SF.txt');
%dataB=load('SF/res_sum_table_banerjee_SF.txt');
dataV=load('random/res_sum_table_exp_vector_random.txt');
i30V = find(dataV(:,1)==30);
i50V = find(dataV(:,1)==50);
i500V = find(dataV(:,1)==500);
i1000V = find(dataV(:,1)==1000);
 
 
figure(22);
hold on;
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30,9),data(i30,8),'k-');
plot(dataV(i30V(1),6),dataV(i30V(1),5),'kv','MarkerSize',12);

plot(data(i50,9),data(i50,8),'g-');
plot(dataV(i50V(1),6),dataV(i50V(1),5),'gv','MarkerSize',12);

plot(data(i500,9),data(i500,8),'b-');
plot(dataV(i500V(1),6),dataV(i500V(1),5),'bv','MarkerSize',12);

plot(data(i1000,9),data(i1000,8),'r-');
plot(dataV(i1000V(1),6),dataV(i1000V(1),5),'rv','MarkerSize',12);

 title('random networks','FontSize', 16);
legend('N=30','Vector(N=30)','N=50','Vector(N=50)','N=500','Vector(N=500)','N=1000','Vector(N=1000)');

figure(3);
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30sp2,6),data(i30sp2,5),'k-.');
hold on;
plot(data(i500sp2,6),data(i500sp2,5),'b--');
plot(data(i1000sp2,6),data(i1000sp2,5),'r-');

figure(4);
xlabel('False Positives','FontSize', 16);
ylabel('True Positives','FontSize', 16);
plot(data(i30sp3,6),data(i30sp3,5),'k-.');
hold on;
plot(data(i500sp3,6),data(i500sp3,5),'b--');bv
plot(data(i1000sp3,6),data(i1000sp3,5),'r-');



