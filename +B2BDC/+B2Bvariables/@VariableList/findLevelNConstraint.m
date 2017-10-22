function varList = findLevelNConstraint(obj,x)
% VARLIST = FINDLEVELNCONSTRAINT(OBJ,X) returns a new
% variableList object with added linear constraints.

%  Created: Sep 26, 2017     Wenyu Li

varList = obj;
[idx,c] = kmeans(x,2);
cc = norm(c(1,:)-c(2,:));
svmModel = fitcsvm(x,idx,'KernelFunction','linear',...
   'Standardize',false);
vec_coef = svmModel.Beta;
id = find(vec_coef,1);
vec_coef = vec_coef/vec_coef(id);
tt = x*vec_coef;
tt_pos = tt(tt>=0);
tt_neg = tt(tt<0);
LB = max(tt_neg);
UB = min(tt_pos);
if UB-LB >= 0.2*cc
   varList= varList.addLinearConstraint(vec_coef',LB,UB);
else
   disp('The outliers are too close to seperate')
end


   


