function alpha=findnegstep(K, W, i, j, A, S, Cp, Cpp, update) 
% This function sets up and solves quadratic equation to find a step length
% for the step C'(i,j)=C'(i,j)-alpha (update=1) or C''(i,j)=C''(i,j)-alpha
% (update  =-1). Notice that it is not necessary to pass A, W and S, but
% only their appropriate elements in the positions (i,j), (j,i), (i,i) and
% (j,j)
if update ==1
  a1=2*(K*W(i,j)-A(i,j))-S(j,i)-S(i,j);
  a2=(W(i,i)*W(j,j)-W(i,j)^2);
  a3=2*W(i,j);
  a4=-2*K*(W(i,i)*W(j,j)+W(i,j)^2);
  a5=2*K*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2);
  maxstep=Cp(i,j);
else
  a1=2*(-K*W(i,j)+A(i,j))-S(i,j)-S(j,i);  
  a2=(-W(i,i)*W(j,j)+W(i,j)^2);
  a3=2*W(i,j);
  a4=2*K*(W(i,i)*W(j,j)+W(i,j)^2);
  a5=2*K*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2);
  maxstep=Cpp(i,j);
end
% If the diaonal element is updated, then the step length needs to be
% halved, since C(i,i) will be updated by alpha twice.
if (i == j)
  maxstep=maxstep/2;
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
  if alpha2<0
    alpha=max(alpha2, -maxstep);
%    display 'possibly nonunique step size'
%    keyboard
  elseif alpha1<0
    alpha=max(alpha1, -maxstep);
  else
   alpha=-maxstep; 
  end
% Make sure that the element C'(i,j) or C''(i,j) remains nonnegative
elseif (-c/b<0)
  alpha=max(-c/b, -maxstep);
else
  alpha=-maxstep;
end
% Check for debuggin purposes. 
% theta=alpha;
% if update ==1 
% fprime=(K*W(i,j)+K*W(j,i)-A(i,j)-A(j,i)-S(i,j)-S(j,i))- ...
%      K*theta/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j))*(-2*(W(i,i)*W(j,j)+W(i,j)^2) ...
%    +2*theta*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2));
% else
%    fprime=(-K*W(i,j)-K*W(j,i)+A(i,j)+A(j,i)-S(i,j)-S(j,i))-K*theta/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)... 
%            +1-2*theta*W(i,j))*(2*(W(i,i)*W(j,j)+W(i,j)^2) ...
%            +2*theta*W(i,j)*(W(i,i)*W(j,j)-W(i,j)^2));
% end
% if fprime > 10^(-6)
%   display 'bad gradient'
%   keyboard
% end
%a=(1-theta*W(i,j))/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
%b=theta*W(j,j)/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
%c=theta*W(i,i)/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
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