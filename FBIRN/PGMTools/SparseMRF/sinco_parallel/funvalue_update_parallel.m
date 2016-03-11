function fchange=funvalue_update_parallel(alpha, K, Ws, i, j, As, Ss, update);

% Modified this function to be more efficient when used with parfor. We avoid sending the whole
% matrices for A,S,W. Instead we replace each with a struct that has params ij,ji,ii,jj. 
% E.g. Ws.ii

% This function computed the update of the objective function if C'(i,j) becomes C'(i,j)+alpha
% (update=1) or C'(,ij)'=C''(i,j)+alpha (update=-1) (alpha can be negative here). 
% Notice that it is not necessary to pass A, W and S, but
% only their appropriate elements in the positions (i,j), (j,i), (i,i) and (j,j)
if update ==1
%  detCnew=detC*(1+alpha*W(i,j))*(1+alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));

 %detratio=(1+alpha*W(i,j))*(1+alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));
  detratio=(1+alpha*Ws.ij)*(1+alpha*Ws.ij-alpha^2*Ws.ii*Ws.jj/(1+alpha*Ws.ij));
 %fchange=K*(log(detratio))-2*alpha*A(i,j)-alpha*(S(i,j)+S(j,i));
  fchange=K*(log(detratio))-2*alpha*As.ij-alpha*(Ss.ij+Ss.ji);
else
%  detCnew=detC*(1-alpha*W(i,j))*(1-alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1+alpha*W(i,j)));

 %detratio=(1-alpha*W(i,j))*(1-alpha*W(i,j)-alpha^2*W(i,i)*W(j,j)/(1-alpha*W(i,j)));
  detratio=(1-alpha*Ws.ij)*(1-alpha*Ws.ij-alpha^2*Ws.ii*Ws.jj/(1-alpha*Ws.ij));
 %fchange=K*(log(detratio))+2*alpha*A(i,j)-alpha*(S(i,j)+S(j,i));
  fchange=K*(log(detratio))+2*alpha*As.ij-alpha*(Ss.ij+Ss.ji);
end
