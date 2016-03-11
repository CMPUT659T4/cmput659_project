function fchange=funvalue_update(alpha, K, W, i, j, A, S, update);
% This function computed the update of the objective function if C'(i,j) becomes C'(i,j)+alpha
% (update=1) or C'(,ij)'=C''(i,j)+alpha (update=-1) (alpha can be negative here). 
% Notice that it is not necessary to pass A, W and S, but
% only their appropriate elements in the positions (i,j), (j,i), (i,i) and (j,j)
if update ==1
%  detCnew=detC*(1+alpha*W(i,j))*(1+alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));
  detratio=(1+alpha*W(i,j))*(1+alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));
  fchange=K*(log(detratio))-2*alpha*A(i,j)-alpha*(S(i,j)+S(j,i));
else
%  detCnew=detC*(1-alpha*W(i,j))*(1-alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));
  detratio=(1-alpha*W(i,j))*(1-alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1-alpha*W(i,j)));
  fchange=K*(log(detratio))+2*alpha*A(i,j)-alpha*(S(i,j)+S(j,i));
end
