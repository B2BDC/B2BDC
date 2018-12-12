function [x,y] = calculatePrediction(ds,QQ,QOIgroup,f1,f2,design,range,Hv)

if nargin < 8
   Hv = [];
end
if isempty(QOIgroup)
   QOIgroup = {1:ds.Length};
end
if ~iscell(QOIgroup)
   QOIgroup = {QOIgroup};
end
nGroup = length(QOIgroup);
if ~iscell(f1)
   f1 = {f1};
end
if length(f1) ~= nGroup
   error('lol');
end
if ~isempty(f1{1})
   nmd = zeros(nGroup,1);
   GroupIndex = zeros(ds.Length,1);
   for j = 1:nGroup
      qoiIndex = QOIgroup{j};
      GroupIndex(qoiIndex) = j;
      nmd(j) = length(f1{j}(design(1,:)));
   end
else
   nmd = [];
   GroupIndex = zeros(ds.Length,1);
end
vars = ds.Variables.Values;
vList = ds.Variables;
n_variable = length(vars);
if isempty(Hv)
   Hv = repmat([-inf inf],n_variable,1);
elseif size(Hv,1) ~= n_variable
   error(' lol ');
end
nMD = sum(nmd);
nPD = 0;
if ~isempty(f2)
   npd = zeros(n_variable,1);
   for j = 1:n_variable
      if ~isempty(f2{j})
         npd(j) = length(f2{j}(design(1,:)));
      end
   end
   nPD = sum(npd);
end
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
LB = [LB; -range*ones(nMD+nPD,1)];
UB = [UB; range*ones(nMD+nPD,1)];
A1 = [];
B1 = [];
Aeq = [];
Beq = [];
units = ds.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
allVarnames = ds.VarNames;
idall = cell(n_units,1);
Qall = cell(n_units,1);
for j = 1:n_units
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   Qall{j} = tmodel.CoefMatrix;
end
tmodel = QQ.model;
[~,~,idpred] = intersect(tmodel.VarNames,allVarnames,'stable');
Qpred = tmodel.CoefMatrix;
predGroup = QQ.group;
spred = QQ.sv;
if predGroup > 0
   predBasis = f1{predGroup}(spred);
end
if nPD > 0
   vBasis = cell(length(idpred),1);
   for j = 1:length(idpred)
      if ~isempty(f2{idpred(j)})
         vBasis{j} = f2{idpred(j)}(spred);
      end
   end
   Apd = [];
   bpd = [];
   for j = 1:n_units
      id1 = idall{j};
      for k = 1:length(id1)
         if ~isempty(f2{id1(k)})
            AA = zeros(2,n_variable+nMD+nPD);
            AA(:,id1(k)) = [1;-1];
            bb = [Hv(id1(k),2); -Hv(id1(k),1)];
            AA(:,n_variable+nMD+sum(npd(1:id1(k)-1))+1:n_variable+nMD+sum(npd(1:id1(k)))) =...
               [f2{id1(k)}(design(j,:)); -f2{id1(k)}(design(j,:))];
            Apd = [Apd; AA];
            bpd = [bpd; bb];
         end
      end
   end
   infFlag = find(isinf(bpd));
   bpd(infFlag) = [];
   Apd(infFlag,:) = [];
   A1 = [A1;Apd];
   B1 = [B1;bpd];
end
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','MaxIter',1500,...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
nStart = 5;
xStart = zeros(2*nStart,n_variable+nMD+nPD);
xStart(:,2:n_variable+1) = vList.makeLHSsample(2*nStart);
% xtest = [];
% load xtest
% if strcmp(CM,'rel')
%    y1 = neq(xtest');
% elseif strcmp(CM,'abs')
%    y1 = neq2(xtest');
% end
% save y1 y1
xmin = cell(nStart,1);
ymin = zeros(nStart,1);
flagmin = zeros(nStart,1);
xmax = cell(nStart,1);
ymax = zeros(nStart,1);
flagmax = zeros(nStart,1);

% xtest =[];
% x = []; xx = [];
% load xtest
% y1.min = funmin(xtest');
% y1.max = funmax(xtest');
% y1.neq = neq(xtest');
% save y1 y1
% load xopt
% load xopt2

for j = 1:nStart
   [xmin{j},ymin(j),flagmin(j)] = fmincon(@funmin,xStart(2*j-1,:)',...
      A1,B1,Aeq,Beq,LB,UB,@neq,opt);
   [xmax{j},ymax(j),flagmax(j)] = fmincon(@funmax,xStart(2*j,:)',...
      A1,B1,Aeq,Beq,LB,UB,@neq,opt);
end
ymin(flagmin<=0) = [];
xmin(flagmin<=0) = [];
ymax(flagmax<=0) = [];
xmax(flagmax<=0) = [];
if ~isempty(ymin)
   [y.min,ID] = min(ymin);
   x.min = xmin{ID};
end
if ~isempty(ymax)
   [y.max,ID] = min(ymax);
   y.max = -y.max;
   x.max = xmax{ID};
end


   function [y,gy] = funmin(x)
      dx = zeros(length(idpred),1);
      for i = 1:length(idpred)
         if ~isempty(vBasis{i})
            dx(i) = vBasis{i}*x(n_variable+nMD+sum(npd(1:idpred(i)-1))+1:n_variable+nMD+sum(npd(1:idpred(i))));
         end
      end
      y = [1;x(idpred)+dx]'*Qpred*[1;x(idpred)+dx]+...
         predBasis*x(n_variable+sum(nmd(1:predGroup-1))+1:n_variable+sum(nmd(1:predGroup)));
      gy = zeros(n_variable+nMD+nPD,1);
      gy(idpred) = 2*(Qpred(2:end,2:end)*(x(idpred)+dx)+Qpred(2:end,1));
      for i = 1:length(idpred)
         if ~isempty(vBasis{i})
            gy(n_variable+nMD+sum(npd(1:idpred(i)-1))+1:n_variable+nMD+sum(npd(1:idpred(i)))) =...
               gy(idpred(i))*vBasis{i};
         end
      end
      gy(n_variable+sum(nmd(1:predGroup-1))+1:n_variable+sum(nmd(1:predGroup))) = predBasis;
   end

   function [y,gy] = funmax(x)
      dx = zeros(length(idpred),1);
      for i = 1:length(idpred)
         if ~isempty(vBasis{i})
            dx(i) = vBasis{i}*x(n_variable+nMD+sum(npd(1:idpred(i)-1))+1:n_variable+nMD+sum(npd(1:idpred(i))));
         end
      end
      y = -[1;x(idpred)+dx]'*Qpred*[1;x(idpred)+dx]-...
         predBasis*x(n_variable+sum(nmd(1:predGroup-1))+1:n_variable+sum(nmd(1:predGroup)));
      gy = zeros(n_variable+nMD+nPD,1);
      gy(idpred) = -2*(Qpred(2:end,2:end)*(x(idpred)+dx)+Qpred(2:end,1));
      for i = 1:length(idpred)
         if ~isempty(vBasis{i})
            gy(n_variable+nMD+sum(npd(1:idpred(i)-1))+1:n_variable+nMD+sum(npd(1:idpred(i)))) =...
               gy(idpred(i))*vBasis{i};
         end
      end
      gy(n_variable+sum(nmd(1:predGroup-1))+1:n_variable+sum(nmd(1:predGroup))) = -predBasis;
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         quadCoef = Qall{i};
         tx = x(id);
         dx = zeros(length(id),1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  dx(ii) = f2{id(ii)}(design(i,:))*...
                     x(n_variable+nMD+sum(npd(1:id(ii)-1))+1:n_variable+nMD+sum(npd(1:id(ii))));
               end
            end
         end
         c(2*i-1) = [1;tx+dx]'*quadCoef*[1;tx+dx]-u;
         c(2*i) =  l-[1;tx+dx]'*quadCoef*[1;tx+dx];
         g(id,2*i-1) = 2*quadCoef(2:end,2:end)*(tx+dx)+2*quadCoef(2:end,1);
         g(id,2*i) = -g(id,2*i-1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  g(n_variable+nMD+sum(npd(1:id(ii)-1))+1:n_variable+nMD+sum(npd(1:id(ii))),2*i-1) =...
                     g(id(ii),2*i-1)*f2{id(ii)}(design(i,:));
                  g(n_variable+nMD+sum(npd(1:id(ii)-1))+1:n_variable+nMD+sum(npd(1:id(ii))),2*i) =...
                     g(id(ii),2*i)*f2{id(ii)}(design(i,:));
               end
            end
         end
         tmpGroup = GroupIndex(i);
         if tmpGroup > 0 
            c(2*i-1) = c(2*i-1) + f1{tmpGroup}(design(i,:))*...
               x(n_variable+sum(nmd(1:tmpGroup-1))+1:n_variable+sum(nmd(1:tmpGroup)));
            c(2*i) = c(2*i) - f1{tmpGroup}(design(i,:))*...
               x(n_variable+sum(nmd(1:tmpGroup-1))+1:n_variable+sum(nmd(1:tmpGroup)));
            g(n_variable+sum(nmd(1:tmpGroup-1))+1:n_variable+sum(nmd(1:tmpGroup)),2*i-1) = f1{tmpGroup}(design(i,:));
            g(n_variable+sum(nmd(1:tmpGroup-1))+1:n_variable+sum(nmd(1:tmpGroup)),2*i) = - f1{tmpGroup}(design(i,:));
         end
      end
      ceq = [];
      geq = [];
   end

end