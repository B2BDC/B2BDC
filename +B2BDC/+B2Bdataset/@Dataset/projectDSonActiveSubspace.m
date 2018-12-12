function newDS = projectDSonActiveSubspace(obj,PCinfo,nPC)

% NEWDS = PROJECTDSONACTIVESUBSPACE(OBJ,PCINFO,NPC) projects the dataset
% onto a lower-dimensional active subspace, specified by PCinfo and nPC.

%  Created: Jan 18, 2018     Wenyu Li

vv = PCinfo.V;
% dd = PCinfo.D;
% [~,id] = sort(dd,'descend');
% vv = vv(:,id);
v1 = vv(:,1:nPC);
x0 = PCinfo.x0;
y0 = x0*v1;
var0 = obj.Variables;
nvar = var0.Length;
H0 = var0.calBound;
newH = zeros(nPC,2);
dH1 = H0(:,2)-x0';
dH2 = x0' - H0(:,1);
for i = 1:nPC
   tv = vv(:,i);
   idPos = tv>=0;
   idNeg = tv<0;
   tt = inf(nvar,1);
   tt(idPos) = dH1(idPos)./tv(idPos);
   tt(idNeg) = -dH2(idNeg)./tv(idNeg);
   newH(i,2) = min(tt);
   tv = -tv;
   idPos = tv>=0;
   idNeg = tv<0;
   tt = inf(nvar,1);
   tt(idPos) = dH1(idPos)./tv(idPos);
   tt(idNeg) = -dH2(idNeg)./tv(idNeg);
   newH(i,1) = -min(tt);
end
varName = cell(nPC,1);
for i = 1:nPC
   varName{i} = ['Active variable ' num2str(i)];
end
newVar = generateVar(varName,newH,zeros(nPC,1));
newDS = generateDataset('Transformed dataset in active subspace');
name0 = {var0.Values.Name}';
obs = obj.calObserve;
bds = obj.calBound;
for i = 1:obj.Length
   unitName = obj.DatasetUnits.Values(i).Name;
   m0 = obj.DatasetUnits.Values(i).SurrogateModel;
   err0 = m0.ErrorStats;
   coef0 = m0.CoefMatrix;
   name1 = m0.VarNames;
   [~,~,id] = intersect(name1,name0,'stable');
   a0 = zeros(nvar,1);
   a0(id) = 2*coef0(2:end,1);
   b0 = zeros(nvar);
   b0(id,id) = coef0(2:end,2:end);
   a1 = 0.5*(a0'+2*x0*b0)*v1;
   b1 = 0.5*v1'*b0*v1;
   coefNew = zeros(nPC+1);
   coefNew(1,2:end) = a1;
   coefNew(2:end,2:end) = b1;
   coefNew = coefNew + coefNew';
   coefNew(1,1) = coef0(1,1)+a0'*x0'+x0*b0*x0';
   mNew = generateModel(coefNew,newVar);
   mNew.ErrorStats = err0;
   unitNew = generateDSunit(unitName,mNew,bds(i,:),obs(i));
   newDS.addDSunit(unitNew);
end
opt = generateOpt('Display',false);
if ~obj.FeasibleFlag
   opt.AddFitError = false;
end
for i = 1:100 
   if newDS.isConsistent(opt)
      return
   else
      newDS.clearConsis;
   end
end
error('The projected dataset is inconsistent');
