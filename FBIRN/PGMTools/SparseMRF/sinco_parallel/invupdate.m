function [Gradp, Gradpp, UpInd, W]=invupdate(theta, Gradp, Gradpp, W, K, i, j, update);
% This function computed the update to W and the gradient matrices when alpha is added to C'(i,j) 
% (update=1) or to C''(i,j) (update=-1). Notice that what is computed is the update 
% UpW, which is a rank-4 matrix obtained from weighted outer products of columns of W.
% It is not necessary to pass A, W, Gradp, Gradpp and S, but
% only  appropriate elements of A, W and S in the positions (i,j), (j,i), (i,i) and
% (j,j) and then pass back a, b and c for the actual matrix update.
if update == 1
 a=-(1+theta*W(i,j))/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
 b=theta*W(j,j)/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
 c=theta*W(i,i)/(theta^2*(W(i,i)*W(j,j)-W(i,j)^2)-1-2*theta*W(i,j));
 UpW=-theta*(a*W(:,i)*W(:,j)'+ c*W(:,j)*W(:,j)'+b*W(:,i)*W(:,i)'+a*W(:,j)*W(:,i)');
else
 a=(1-theta*W(i,j))/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
 b=theta*W(j,j)/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
 c=theta*W(i,i)/(theta^2*(-W(i,i)*W(j,j)+W(i,j)^2)+1-2*theta*W(i,j));
 UpW=theta*(a*W(:,i)*W(:,j)'+ c*W(:,j)*W(:,j)'+b*W(:,i)*W(:,i)'+a*W(:,j)*W(:,i)');
end
W=W+UpW;
Gradp=Gradp+K*UpW;
Gradpp=Gradpp-K*UpW;
% UpInd should be nonzero wherever matrix W needed to be updated.
UpInd=(abs(UpW)>0);