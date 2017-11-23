function ds = generateSubDS(obj,vv,nPC,xAve,sVar)

% generate a dataset in the truncated principal direction space.

%  Created: June 7, 2017     Wenyu Li


ds = generateDataset('Trunced dataset');
units = obj.DatasetUnits.Values;
n_unit = length(units);
allName = obj.VarNames;
v1 = vv(:,1:nPC);
% v2 = vv(:,nPC+1:end);
if size(xAve,1) == 1
   xAve = xAve';
end
% yAve = v2*xAve;
% y2 = v2'*xAve;
bds = obj.calBound;
obs = obj.calObserve;
for i = 1:n_unit
   tModel = units(i).SurrogateModel;
   tName = tModel.VarNames;
   [~,~,id] = intersect(allName,tName,'stable');
   Coef_old = tModel.CoefMatrix;
   Coef_old = Coef_old([1;id+1],[1;id+1]);
   Coef_new = zeros(nPC+1);
   a = 2*Coef_old(2:end,1);
   Q = Coef_old(2:end,2:end);
   Coef_new(1,1) = Coef_old(1,1)+a'*xAve+xAve'*Q*xAve;
   a_new = a'*v1 + 2*xAve'*Q*v1;
   Coef_new(1,2:end) = 0.5*a_new;
   Coef_new(2:end,1) = 0.5*a_new';
   Coef_new(2:end,2:end) = v1'*Q*v1;
   Coef_new = 0.5*(Coef_new+Coef_new');
   new_model = generateModel(Coef_new,sVar);
   newDSunit = generateDSunit(units(i).Name,new_model,bds(i,:),obs(i));
   ds.addDSunit(newDSunit);
end


