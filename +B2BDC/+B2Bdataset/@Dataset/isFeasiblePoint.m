function [isFeas,xFea] = isFeasiblePoint(ds, point)
%   [ISFEAS,xFea] = ISFEASIBLEPOINT(DS, POINT) returns a logical true if POINT is 
%   feasible for the dataset DS. POINT is a nPoint-by-nVariables vector of 
%   the points which you want to examine feasibility. If any points are
%   infeasible, a logical false is returned.

tolY = 0;
nVar = ds.Variables.Length;
% if size(point,2) ~= nVar
%    point = point';
% end
% if size(point,2) ~= nVar
%    error('Wrong input point size')
% end
n_units = ds.Length;
nn = nargout;
[nPoint,nD] = size(point);
H = ds.Variables.calBound;
if ds.ModelDiscrepancyFlag
   nMD = ds.ModelDiscrepancy.Variables.Length;
   xMD = ds.ModelDiscrepancy.FeasiblePoint;
   HMD = ds.ModelDiscrepancy.Variables.calBound;
   if isempty(xMD)
      xMD = 0.5*(HMD(:,1)+HMD(:,2));
   end
   if nD == nVar
      H = [H; HMD];
   end
else
   xMD = [];
   nMD = 0;
end
if ds.ParameterDiscrepancyFlag
   nPD = ds.ParameterDiscrepancy.Variables.Length;
   xPD = ds.ParameterDiscrepancy.FeasiblePoint;
   HPD = ds.ParameterDiscrepancy.Variables.calBound;
   if isempty(xPD)
      xPD = 0.5*(HPD(:,1)+HPD(:,2));
   end
   if nD == nVar
      H = [H; HPD];
   end
else
   xPD = [];
   nPD = 0;
end
if nD ~= nVar && nD ~= nVar+nMD+nPD
   point = point';
end
if nD == nVar+nMD+nPD
   xFea = point;
elseif nD == nVar
   xFea = [];
else
   error('Invalid input sample point size');  
end
isFeas = false(nPoint,1);
dsBounds = ds.calBound;
dbds = diff(dsBounds,[],2);
if ds.FeasibleFlag
   for j = 1:ds.Length
      abE = ds.DatasetUnits.Values(j).SurrogateModel.ErrorStats.absMax;
      dsBounds(j,1) = dsBounds(j,1) - abE;
      dsBounds(j,2) = dsBounds(j,2) + abE;
   end
end
dsBounds(:,1) = dsBounds(:,1) - 0.5*dbds*tolY;
dsBounds(:,2) = dsBounds(:,2) + 0.5*dbds*tolY;
% ic = 0;
nmax = 10^5;
n1 = floor(nPoint/nmax);
n2 = mod(nPoint,nmax);
xid = 1:nmax;
[idall,Qall,~,~,APD,bPD] = ds.getQ_RQ_expansion;
ACM = [zeros(size(APD,1),1) APD];
for i1 = 1:n1
   b1 = ds.Variables.isFeasiblePoint(point(xid,1:nVar));
   ID = find(b1);
   if nD == nVar+nMD+nPD
      dsEval = ds.eval(point(xid(b1),:));
      a1 = dsEval > repmat(dsBounds(:,1)',sum(b1),1);
      a2 = dsEval < repmat(dsBounds(:,2)',sum(b1),1);
      b2 = all(a1');
      b3 = all(a2');
      if ~isempty(APD)
         a3 = APD*point(xid(b1),:)' < repmat(bPD,1,sum(b1));
         b4 = all(a3);
         b5 = b2 & b3 & b4;
      else
         b5 = b2 & b3;
      end
      ID(~b5) = [];
   else
      [newpoint, pid] = findFeasibleDiscrepancy(point(xid(b1),:));
      if nn > 1
         xFea = [xFea; newpoint(pid,:)];
      end
      ID(~pid) = [];
   end
   isFeas(ID+(i1-1)*nmax) = true;
   point(1:nmax,:) = [];
end
if n2 > 0
   xid = 1:n2;
   b1 = ds.Variables.isFeasiblePoint(point(xid,1:nVar));
   ID = find(b1); 
   if nD == nVar+nMD+nPD
      dsEval = ds.eval(point(xid(b1),:));
      a1 = dsEval > repmat(dsBounds(:,1)',sum(b1),1);
      a2 = dsEval < repmat(dsBounds(:,2)',sum(b1),1);
      b2 = all(a1');
      b3 = all(a2');
      if ~isempty(APD)
         a3 = APD*point(xid(b1),:)' < repmat(bPD,1,sum(b1));
         b4 = all(a3);
         b5 = b2 & b3 & b4;
      else
         b5 = b2 & b3;
      end
      ID(~b5) = []; 
   else
      [newpoint, pid] = findFeasibleDiscrepancy(point(xid(b1),:));
      if nn > 1
         xFea = [xFea; newpoint(pid,:)];
      end
      ID(~pid) = [];
   end
   isFeas(ID+n1*nmax) = true;
end
if nD == nVar+nMD+nPD
   xFea = xFea(isFeas,:);
end

   function [newx,flag] = findFeasibleDiscrepancy(x)
      nx = size(x,1);
      flag = true(nx,1);
      xStart = [zeros(nx,1) x repmat([xMD ;xPD]',nx,1)];
      newx = xStart(:,2:end);
      opt = optimoptions('fmincon','Display','none','GradObj','on',...
         'GradConstr','on','Algorithm','interior-point','MaxIter',3000,...
         'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
      Aeq = [zeros(nVar,1) eye(nVar) zeros(nVar,nMD+nPD)];
      for ii = 1:nx
         [x0,y0,exitflag] = fmincon(@funxmin,xStart(ii,:)',...
            ACM,bPD,Aeq,x(ii,:)',[-Inf;H(:,1)],[Inf;H(:,2)],@neq,opt);
         if exitflag >= 0 && y0 < 0
            newx(ii,nVar+1:end) = x0(nVar+2:end);
         else
            flag(ii) = false;
         end
      end
   end
   
   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1 ; zeros(nVar+nMD+nPD,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(nVar+nMD+nPD+1,2*n_units);
      for i = 1:n_units
         id = idall{i};
         quadCoef = Qall{i};
         c(2*i-1) = [1;x(id+1)]'*quadCoef*[1;x(id+1)];
         g(id+1,2*i-1) = 2*(quadCoef(2:end,2:end)*x(id+1)+quadCoef(2:end,1));
      end
      c(2:2:end) = -c(1:2:end);
      c(1:2:end) = c(1:2:end)+x(1)-dsBounds(:,2);
      c(2:2:end) = c(2:2:end)+x(1)+dsBounds(:,1);
      g(:,2:2:end) = -g(:,1:2:end);
      g(1,:) = 1;
      ceq = [];
      geq = [];
   end
end


