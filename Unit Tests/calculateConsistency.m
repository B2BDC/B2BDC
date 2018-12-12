function [x,y,lambda] = calculateConsistency(ds,QOIgroup,f1,f2,design,range,Hv,CM)

if nargin < 8
   CM = 'rel';
end
if nargin < 7
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
if nPD > 0
   Apd = [];
   bpd = [];
   for j = 1:n_units
      id1 = idall{j};
      for k = 1:length(id1)
         if ~isempty(f2{id1(k)})
            AA = zeros(2,n_variable+nMD+nPD+1);
            AA(:,id1(k)+1) = [1;-1];
            bb = [Hv(id1(k),2); -Hv(id1(k),1)];
            AA(:,n_variable+1+nMD+sum(npd(1:id1(k)-1))+1:n_variable+1+nMD+sum(npd(1:id1(k)))) =...
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
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','on');
nStart = 5;
xStart = zeros(nStart,n_variable+nMD+nPD+1);
xStart(:,2:n_variable+1) = vList.makeLHSsample(nStart);
% xtest = [];
% load xtest
% if strcmp(CM,'rel')
%    y1 = neq(xtest');
% elseif strcmp(CM,'abs')
%    y1 = neq2(xtest');
% end
% save y1 y1
for j = 1:nStart
   if strcmp(CM,'rel')
      [x0,yin_result,~,~,ilam] = fmincon(@funxmin,xStart(j,:)',...
         A1,B1,Aeq,Beq,[-Inf;LB],[Inf;UB],@neq,opt);
   elseif strcmp(CM,'abs')
      [x0,yin_result,~,~,ilam] = fmincon(@funxmin,xStart(j,:)',...
         A1,B1,Aeq,Beq,[-Inf;LB],[Inf;UB],@neq2,opt);
   end
   if yin_result<0
      break
   end
end
x = x0(2:end);
y = -yin_result;
lambda = ilam;




   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1;zeros(n_variable+nMD+nPD,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD+1,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         quadCoef = Qall{i};
         tx = x(id+1);
         dx = zeros(length(id),1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  dx(ii) = f2{id(ii)}(design(i,:))*...
                     x(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))));
               end
            end
         end
         c(2*i-1) = [1;tx+dx]'*quadCoef*[1;tx+dx]-(1-x(1))*u-x(1)*d;
         c(2*i) =  (1-x(1))*l+x(1)*d-[1;tx+dx]'*quadCoef*[1;tx+dx];
         g(id+1,2*i-1) = 2*quadCoef(2:end,2:end)*(tx+dx)+2*quadCoef(2:end,1);
         g(id+1,2*i) = -g(id+1,2*i-1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  g(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))),2*i-1) =...
                     g(id(ii)+1,2*i-1)*f2{id(ii)}(design(i,:));
                  g(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))),2*i) =...
                     g(id(ii)+1,2*i)*f2{id(ii)}(design(i,:));
               end
            end
         end
         tmpGroup = GroupIndex(i);
         if tmpGroup > 0 
            c(2*i-1) = c(2*i-1) + f1{tmpGroup}(design(i,:))*...
               x(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)));
            c(2*i) = c(2*i) - f1{tmpGroup}(design(i,:))*...
               x(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)));
            g(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)),2*i-1) = f1{tmpGroup}(design(i,:));
            g(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)),2*i) = - f1{tmpGroup}(design(i,:));
         end
         g(1,2*i-1) = u-d;
         g(1,2*i) = d-l; 
      end
      ceq = [];
      geq = [];
   end

   function [c,ceq,g,geq] = neq2(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD+1,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         quadCoef = Qall{i};
         tx = x(id+1);
         dx = zeros(length(id),1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  dx(ii) = f2{id(ii)}(design(i,:))*...
                     x(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))));
               end
            end
         end
         c(2*i-1) = [1;tx+dx]'*quadCoef*[1;tx+dx]-u+x(1);
         c(2*i) =  l+x(1)-[1;tx+dx]'*quadCoef*[1;tx+dx];
         g(id+1,2*i-1) = 2*quadCoef(2:end,2:end)*(tx+dx)+2*quadCoef(2:end,1);
         g(id+1,2*i) = -g(id+1,2*i-1);
         if nPD > 0
            for ii = 1:length(id)
               if ~isempty(f2{id(ii)})
                  g(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))),2*i-1) =...
                     g(id(ii)+1,2*i-1)*f2{id(ii)}(design(i,:));
                  g(n_variable+nMD+1+sum(npd(1:id(ii)-1))+1:n_variable+nMD+1+sum(npd(1:id(ii))),2*i) =...
                     g(id(ii)+1,2*i)*f2{id(ii)}(design(i,:));
               end
            end
         end
         tmpGroup = GroupIndex(i);
         if tmpGroup > 0
            c(2*i-1) = c(2*i-1) + f1{tmpGroup}(design(i,:))*...
               x(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)));
            c(2*i) = c(2*i) - f1{tmpGroup}(design(i,:))*...
               x(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)));
            g(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)),2*i-1) = f1{tmpGroup}(design(i,:));
            g(n_variable+1+sum(nmd(1:tmpGroup-1))+1:n_variable+1+sum(nmd(1:tmpGroup)),2*i) = - f1{tmpGroup}(design(i,:));
         end
         g(1,2*i-1) = 1;
         g(1,2*i) = 1;
      end
      ceq = [];
      geq = [];
   end
end