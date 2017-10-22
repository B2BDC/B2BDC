function [varList,constraint] = findLevel2Constraint(obj,x,nc)
% VARLIST = FINDLEVEL2CONSTRAINT(OBJ,X,NC) returns a new
% variableList object with added linear constraints. The number of
% constraints is specified by nc and the linear constraints only involve 2
% variables.

%  Created: Sep 25, 2017     Wenyu Li

varList = obj;
[~,nVar] = size(x);
[x_cov,xlabel,cc] = findCorrelation(x);
[~,id] = sort(x_cov(:),'descend');
A = zeros(nc,nVar);
LB = zeros(nc,1);
UB = zeros(nc,1);
it = 1;
count = 1;
while it <= nc
   [ii,jj] = ind2sub(size(x_cov),id(count));
   tmpx = x(:,[ii,jj]);
   svmModel = fitcsvm(tmpx,xlabel{ii,jj},'KernelFunction','linear',...
    'Standardize',false);
   vec_coef = svmModel.Beta;
   vec_coef = vec_coef/vec_coef(1);
   tt = tmpx*vec_coef;
   tt_pos = tt(tt>=0);
   tt_neg = tt(tt<0);
   if isempty(tt_pos) || isempty(tt_neg)
      count = count+1;
      continue
   end
   lb = max(tt_neg);
   ub = min(tt_pos);
   if ub-lb >= 0.2*cc(ii,jj)
      count = count+1;
      continue
   else
      LB(it) = lb;
      UB(it) = ub;
      A(it,[ii,jj]) = vec_coef;
      count = count+1;
      it = it+1;
   end
end
id = any(A,2);
A = A(id,:);
LB = LB(id);
UB = UB(id);
varList= varList.addLinearConstraint(A,LB,UB);
constraint.A = A;
constraint.LB = LB;
constraint.UB = UB;

   function [score,xlabel,cc] = findCorrelation(x)
      score = zeros(nVar);
      xlabel = cell(nVar,nVar);
      cc = zeros(nVar);
      for i1 = 1:nVar-1
         for i2 = i1+1:nVar 
            tmpx = x(:,[i1,i2]);
            [idx,c] = kmeans(tmpx,2);
            cc(i1,i2) = norm(c(1,:)-c(2,:));
            xlabel{i1,i2} = idx;
            x1 = tmpx(idx==1,:);
            x2 = tmpx(idx==2,:);
            n1 = size(x1,1);
            n2 = size(x2,1);
            score(i1,i2) = n1*corr(x1(:,1),x1(:,2))+n2*corr(x2(:,1),x2(:,2));
         end
      end
      score = abs(score)/(n1+n2);
   end
   
end

   


