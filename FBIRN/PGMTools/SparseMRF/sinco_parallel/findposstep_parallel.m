%function alpha=findposstep(K, W, i, j, A, S, update)  
function alpha=findposstep_parallel(K, Ws, i, j, As, Ss, update)    

% Modified this function to be more efficient when used with parfor. We avoid sending the whole
% matrices for A,S,W. Instead we replace each with a struct that has params ij,ji,ii,jj. E.g. Ws.ii

% This function sets up and solves quadratic equation to find a step length
% for the step C'(i,j)=C'(i,j)+alpha (update=1) or C''(i,j)=C''(i,j)+alpha
% (update  =-1). Notice that it is not necessary to pass A, W and S, but
% only their appropriate elements in the positions (i,j), (j,i), (i,i) and (j,j)
if update ==1
 %a1=2*(K*W(i,j)-A(i,j))-S(j,i)-S(i,j);
  a1=2*(K*Ws.ij-As.ij)-Ss.ji-Ss.ij;
 %a2=(W(i,i)*W(j,j)-W(i,j)^2);
  a2=(Ws.ii*Ws.jj-Ws.ij^2);
 %a3=2*W(i,j);
  a3=2*Ws.ij;
 %a4=-2*K*(W(i,i)*W(j,j)+W(i,j)^2);
  a4=-2*K*(Ws.ii*Ws.jj+Ws.ij^2);
 %a5=2*K*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2);
  a5=2*K*Ws.ij*(Ws.ii*Ws.jj-Ws.ij^2);
else
 %a1=2*(-K*W(i,j)+A(i,j))-S(i,j)-S(j,i);  
  a1=2*(-K*Ws.ij+As.ij)-Ss.ij-Ss.ji;  
 %a2=(-W(i,i)*W(j,j)+W(i,j)^2);
  a2=(-Ws.ii*Ws.jj+Ws.ij^2);
 %a3=2*W(i,j);
  a3=2*Ws.ij;
 %a4=2*K*(W(i,i)*W(j,j)+W(i,j)^2);
  a4=2*K*(Ws.ii*Ws.jj+Ws.ij^2);
 %a5=2*K*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2);
  a5=2*K*Ws.ij*(Ws.ii*Ws.jj-Ws.ij^2);
end
a=a2*a1-a5;
b=-a3*a1-a4;
c=-update*a1;
if (abs(a)>10^(-10)) 
  D=b^2-4*a*c;
  if D<0 
    display 'negative discriminant'
    keyboard
  end
  alpha1=min((-b-sqrt(D))/(2*a),(-b+sqrt(D))/(2*a));
  alpha2=max((-b-sqrt(D))/(2*a),(-b+sqrt(D))/(2*a));
  if alpha1>=0
    alpha=alpha1;
%    display 'possibly an unbounded direction!!!'
%    keyboard
  elseif alpha2>=0
    alpha=alpha2;
  else
    display 'unbounded direction!!!'
    keyboard
  end
elseif (-c/b>0)
  alpha=-c/b;
else
  display 'unbounded direction!!!'
  keyboard
end
% check derivative if debgging is necessary
% theta=alpha;
%keyboard
%a=-(1+theta*W(i,j))/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
%b=theta*W(j,j)/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
%c=theta*W(i,i)/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
%Wnew=W-theta*(a*W(:,i)*W(:,j)'+ c*W(:,j)*W(:,j)'+b*W(:,i)*W(:,i)'+a*W(:,j)*W(:,i)');
%[p,p]=size(W);
%E=eye(p,p);
%Cnew=Cp+theta*E(:,i)*E(j,:)+theta*E(:,j)*E(i,:);
%fprime-sum(sum((I1*I2).*Wnew))
%Wnew-inv(Cnew)
%f=K*log(det(Cp-Cpp))-sum(sum(A.*(Cp-Cpp)))-sum(sum(S.*(Cp+Cpp)));
%fnew=K*log(det(Cnew-Cpp))-sum(sum(A.*(Cnew-Cpp)))-sum(sum(S.*(Cnew+Cpp)));
%(fnew-f)/theta
%fprime
%keyboard
