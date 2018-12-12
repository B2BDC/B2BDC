function varList = findConstraint_quadratic_old(obj,xg,xb,misRate)
% VARLIST = FINDLEVEL2CONSTRAINT(OBJ,XG,XB,MISRATE) returns a new variableList
% object with added quadratic constraints. The number of principal
% directions used is specified by the misclassification rate misRate

%  Created: Oct 13, 2017     Wenyu Li

[n1,nVar] = size(xg);
% n2 = size(xb,1);
if nargin < 4 || misRate < 0
   misRate = 0;
end
m1 = 1;
cuttingDegree = 0.1;
npc = 10;
np = nVar;
% xc = xg - repmat(mean(xg),n1,1);
% x1 = [ones(n1,1) xg];
% x2 = [ones(n2,1) xb];
[V,D] = eig(xg' * xg/n1);
[~,idpca] = sort(diag(D),'descend');
V = V(:,idpca);
% yc = mean(xg)*V;
yg = xg*V ;
% yb = xb*V ;
% sigma = max(max(abs(yg)),quantile(abs(yb),cuttingDegree));
sigma = max(abs(yg));
% sigma = 3*sqrt(diag(D));
D = diag(1./(sigma.^2));
% ypc = meam(xg)*V;
% X = findEllipse;
while m1>misRate && np > 1
   Q = findEllipsePCA(np);
%    UB = 1;
   tt = diag(xb*Q*xb');
   t_shift = quantile(tt,cuttingDegree);
   UB = max(t_shift,1);
   t2 = diag(xg*Q*xg');
   flag = find(t2<=UB);
   m1 = (n1-length(flag))/n1;
   if m1 <= misRate
      break
   end
   if np <= npc
      np = np-1;
   else
      np = max(round(0.5*np),npc);
   end
end
epslon = 1e-6;
Q = Q+epslon*eye(nVar);
id = randperm(length(flag),1);
x0 = xg(flag(id),:);
varList = obj.addQuadraticConstraint(Q,UB,x0);

      
   function X = findEllipse(np)
      vv = V(:,1:np);
      K.l = n1;
      K.s = np;
      A = zeros(n1,n1+np^2);
      for i = 1:n1
         tz = vv'*xg(i,:)'*xg(i,:)*vv;
         A(i,n1+1:end) = tz(:);
      end
      A(:,1:n1) = eye(n1);
      b = ones(n1,1);
      c = zeros(n1+np^2,1);
      tn = eye(np);
      c(n1+1:end) = tn(:);
      pars.fid = 0;
      tic;
      [xopt,~,info] = sedumi(A,b,c,K,pars);
      toc
      X = mat(xopt(n1+1:end));
      X = vv*X*vv';
      X = 0.5*(X+X');
   end

   function X = findEllipsePCA(np)
      VV = V(:,1:np);
      DD = D(1:np,1:np);
      X = VV*DD*VV';
      X = 0.5*(X+X');
   end



end