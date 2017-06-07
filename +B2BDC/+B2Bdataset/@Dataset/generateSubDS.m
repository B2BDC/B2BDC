function ds = generateSubDS(obj,vv,nPC,xAve,sVar)

% generate a dataset in the truncated principal direction space.

%  Created: June 7, 2017     Wenyu Li


ds = generateDataset('Trunced dataset');
units = obj.DatasetUnits.Values;
n_unit = length(units);
allName = obj.VarNames;
v1 = vv(:,1:nPC);
v2 = vv(:,nPC+1:end);
if size(xAve,1) == 1
   xAve = xAve';
end
% yAve = v2*xAve;
y2 = v2'*xAve;
bds = obj.calBound;
obs = obj.calObserve;
for i = 1:n_unit
   tModel = units(i).SurrogateModel;
   tName = tModel.VarNames;
   [~,~,id] = intersect(tName,allName,'stable');
   Coef_old = tModel.CoefMatrix;
   Coef_new = zeros(nPC+1);
   Coef_new(1,1) = Coef_old(1,1);
   c_old = Coef_old(1,2:end);
   c_new = c_old*v1(id,:);
   Coef_new(1,1) = Coef_new(1,1) + c_old*v2(id,:)*y2;
   q_old = Coef_old(2:end,2:end);
   q_new = v1(id,:)'*q_old*v1(id,:);
   Coef_new(1,1) = Coef_new(1,1) + y2' * v2(id,:)' * q_old * v2(id,:) * y2;
   Coef_new(1,2:end) = c_new;
   Coef_new(2:end,1) = c_new';
   Coef_new(2:end,2:end) = q_new;
   Coef_new = 0.5*(Coef_new+Coef_new');
   new_model = generateModel(Coef_new,sVar);
   newDSunit = generateDSunit(units(i).Name,new_model,bds(i,:),obs(i));
   ds.addDSunit(newDSunit);
end
