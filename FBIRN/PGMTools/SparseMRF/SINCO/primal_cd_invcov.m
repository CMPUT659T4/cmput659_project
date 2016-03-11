function [Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, A, S, p, K);
precision = 10^(-7);
f=fstart;
W=Wstart;
tol=10^(-6);
Cp=Cstart.*(Cstart>0);
Cpp=-Cstart.*(Cstart<0);
Gradp=-S+K*W-A;
Gradpp=-S-K*W+A;
C=[vec(Cp);vec(Cpp)];
G=[vec(Gradp);vec(Gradpp)];
I1=find(C>10^(-6));
I2=find(C<=10^(-6));
UpInd=ones(p,p);
fchange=zeros(p,p);
alphas=zeros(p,p);
updates=zeros(p,p);
while (norm(G(I1))>tol || max(G(I2))>tol)
%find the largest positive derivative
  alphamax=0;
  fmax=f;
  for i=1:p
      for j=1:i
        if (UpInd(i,j) || UpInd(i,i) || UpInd(j,j))
          alphas(i,j)=0;
          if (Gradp(i,j)+Gradp(j,i) > tol) & (Cpp(i,j)<=tol)
            updates(i,j)=1;
            alphas(i,j)=findposstep(K, W, i, j, A, S, updates(i,j));    
          elseif (Gradpp(i,j)+Gradpp(j,i) > tol) & (Cp(i,j)<=tol)
            updates(i,j)=-1;
            alphas(i,j)=findposstep(K, W, i, j, A, S, updates(i,j)); 
          elseif (Gradp(i,j)+Gradp(j,i)<-tol) & (Cp(i,j)> tol )
            updates(i,j)=1;
            alphas(i,j)=findnegstep(K, W, i, j, A, S, Cp, Cpp, updates(i,j)); 
          elseif (Gradpp(i,j)+Gradpp(j,i)<-tol) & (Cpp(i,j)> tol) 
            updates(i,j)=-1;
            alphas(i,j)=findnegstep(K, W, i, j, A, S, Cp, Cpp, updates(i,j)); 
          end
          if ( abs(alphas(i,j))> 10^(-6))
%         Cptemp=Cp;
%         Cpptemp=Cpp;
%         if update==1
%           Cptemp(i,j)=Cptemp(i,j)+alpha;
%           Cptemp(j,i)=Cptemp(j,i)+alpha;
%         else
%           Cpptemp(i,j)=Cpptemp(i,j)+alpha;
%           Cpptemp(j,i)=Cpptemp(j,i)+alpha;
%         end
           fchange(i,j)=funvalue_update(alphas(i,j), f, K,W, i, j, A, S, updates(i,j));
%         fnewtest=K*log(det(Cptemp-Cpptemp))-sum(sum(A.*(Cptemp-Cpptemp)))-sum(sum(S.*(Cptemp+Cpptemp)));
%         if (abs(fnew-fnewtest)/(1+abs(fnew)) >10^(-6))
%             display 'bad function update'
%            keyboard
%         end
          else
           fchange(i,j)=0;
          end
        end
        fnew=f+fchange(i,j);
        if fnew > fmax
          fmax=fnew;
          imax=i;
          jmax=j;
          alphamax=alphas(i,j);
          updatemax=updates(i,j);
        end
      end
  end
  if alphamax==0 || max(max(fchange))/(abs(f)+1)<precision;
     display 'no more steps' 
%     keyboard
     break
  end
  update=updatemax;
  alpha=alphamax;
  i=imax;
  j=jmax;
  if update==1
   Cp(i,j)=Cp(i,j)+alpha;
   Cp(j,i)=Cp(j,i)+alpha;
  else
   Cpp(i,j)=Cpp(i,j)+alpha;
   Cpp(j,i)=Cpp(j,i)+alpha;
  end
  Wold=W;
  Gradpold=Gradp;
  Gradppold=Gradpp;
  f=fmax;
  %f
%  aaa=f-(K*log(det(Cp-Cpp))-sum(sum(A.*(Cp-Cpp)))-sum(sum(S.*(Cp+Cpp))));
%  if ( abs(aaa) > 10^(-8))
%    keyboard
%  end
%  f=funvalue_update(alpha, f, K,W, i, j, A, S, update)
  [Gradp, Gradpp, UpInd, W]=invupdate(alpha, Gradp, Gradpp, W, K, i, j, update);
%  Gradp=K*W-A-S;
%  Gradpp=-K*W+A-S;
  C=[vec(Cp);vec(Cpp)];
  G=[vec(Gradp);vec(Gradpp)];
  I1=find(C>10^(-6));
  I2=find(C<=10^(-6));
%  f=K*log(det(Cp-Cpp))-sum(sum(A.*(Cp-Cpp)))-sum(sum(S.*(Cp+Cpp)));
end
%fval(iter)=K*log(det(Cp-Cpp))-sum(sum(A.*(Cp-Cpp)))-sum(sum(S.*(Cp+Cpp)))
Csol=Cp-Cpp;
Wsol=W;
fsol=f;

if fsol< fstart
       display 'negative improvement'
       keyboard
end