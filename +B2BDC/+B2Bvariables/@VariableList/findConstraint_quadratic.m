function varList = findConstraint_quadratic(obj,xg,xb,vs)
% VARLIST = FINDLEVEL2CONSTRAINT(OBJ,XG,XB,MISRATE) returns a new variableList
% object with added quadratic constraints. The number of principal
% directions used is specified by the misclassification rate misRate

%  Created: Oct 13, 2017     Wenyu Li

[n1,~] = size(xg);
% n2 = size(xb,1);
if nargin < 4 || vs < 0
   vs = 0;
end
misRate = 0.05;
xc = xg - repmat(mean(xg),n1,1);
x1 = [ones(n1,1) xg];
% x2 = [ones(n2,1) xb];
[V,D] = eig(xc' * xc/n1);
[dd,idpca] = sort(diag(D),'descend');
V = V(:,idpca);
npc = sum(dd>=vs*dd(1));
yc = mean(xg)*V;
yg = xg*V - yc;
% yb = xb*V - yc;
% sigma = max(max(abs(yg)),quantile(abs(yb),cuttingDegree));
sigma = max(abs(yg));
% sigma = 3*sqrt(diag(D));
D = diag(1./(sigma.^2));
% X = findEllipse;
Q = findEllipsePCA(npc);
%    UB = 1;
tt = diag(x1*Q*x1');
t_shift = quantile(tt,1-misRate);
UB = max(t_shift,1);
[~,id] = min(tt);
% t2 = diag(x1*Q*x1');
% epslon = 1e-6;
% Q(2:end,2:end) = Q(2:end,2:end)+epslon*eye(nVar);
% id = randperm(length(flag),1);
x0 = xg(id,:);
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
      yyc = yc(1:np);
      VV = V(:,1:np);
      DD = D(1:np,1:np);
      X = [yyc*DD*yyc', -yyc*DD*VV'; -VV*DD*yyc', VV*DD*VV'];
      X = 0.5*(X+X');
   end



end