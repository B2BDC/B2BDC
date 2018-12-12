function [QOIrange, QOISensitivity, xOpt] = preQOInn(obj,QOIobj,B2Bopt,rflag)
% subfunction to calculate QOI inner bound with fmincon

%  Created: Oct 10, 2016     Wenyu Li

vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
xLB = [vars.LowerBound]';
xUB = [vars.UpperBound]';
dx = xUB-xLB;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   n1 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub-Alb)./sum(abs(A0),2);
   dA = Aub - Alb;
   ieq = (eqTest <= tolerance);
   if any(ieq)
      A1 = [A0(~ieq,:); -A0(~ieq,:)];
      B1 = [Aub(~ieq); -Alb(~ieq)];
      Aeq = A0(ieq,:);
      Beq = 0.5*(Alb(ieq)+Aub(ieq));
   else
      A1 = [A0; -A0];
      Aeq = [];
      B1 = [Aub; -Alb];
      Beq = [];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
allVarnames = obj.VarNames;
[~,~,id0] = intersect(QOIobj.VarNames,allVarnames,'stable');
fg0 = QOIobj.NetGradientFcn;
units = obj.DatasetUnits.Values;
if rflag
   n_units = length(units)-1;
else
   n_units = length(units);
end
abE = zeros(length(units),1);
bds = obj.calBound;
bds = bds(1:n_units,:);
if B2Bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
UB = bds(:,2)+abE;
LB = bds(:,1)-abE;
bds = diff(bds,[],2);
idall = cell(n_units,1);
fg = cell(n_units,1);
for j = 1:n_units
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   fg{j} = tmodel.NetGradientFcn;
end
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','MaxIter',5000,...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
if B2Bopt.Display
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
nStart = B2Bopt.LocalStart;
xStart = vList.makeUniformsample(2*nStart);
xStart(1:2,:) = repmat((obj.FeasiblePoint)',2,1);
infomin = cell(nStart,2);
infomax = cell(nStart,2);
Qmin = zeros(nStart,1);
Qmax = zeros(nStart,1);
flagmin = zeros(nStart,1);
flagmax = zeros(nStart,1);
for j = 1:nStart
   [infomin{j,1},Qmin(j),flagmin(j),~,infomin{j,2}] = fmincon(@funxmin,xStart(2*j-1,:)',A1,B1,Aeq,...
      Beq,xLB,xUB,@neq,opt);
   [infomax{j,1},Qmax(j),flagmax(j),~,infomax{j,2}] = fmincon(@funxmax,xStart(2*j,:)',A1,B1,Aeq,...
      Beq,xLB,xUB,@neq,opt);
end
infomin(flagmin<=0,:) = [];
infomax(flagmax<=0,:) = [];
Qmin(flagmin<=0) = [];
Qmax(flagmax<=0) = [];
[Qmin,idmin] = min(Qmin);
[Qmax,idmax] = min(Qmax);
Qmax = -Qmax;
xmin = infomin{idmin,1};
ilam_min = infomin{idmin,2};
xmax = infomax{idmax,1};
ilam_max = infomax{idmax,2};
s.min.expu = -ilam_min.ineqnonlin(1:2:2*n_units).*bds;
s.min.expl = -ilam_min.ineqnonlin(2:2:2*n_units).*bds;
s.min.varu = -ilam_min.upper.*dx;
s.min.varl = -ilam_min.lower.*dx;
if ~isempty(A1)
   s.min.linu = zeros(n1,1);
   s.min.linl = zeros(n1,1);
   s.min.linu(~ieq) = -ilam_min.ineqlin(1:n1).*dA(~ieq);
   s.min.linl(~ieq) = -ilam_min.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.min.linu = [];
   s.min.linl = [];
end

s.max.expu = ilam_max.ineqnonlin(1:2:2*n_units).*bds;
s.max.expl = ilam_max.ineqnonlin(2:2:2*n_units).*bds;
s.max.varu = ilam_max.upper.*dx;
s.max.varl = ilam_max.lower.*dx;
if ~isempty(A1)
   s.max.linu = zeros(n1,1);
   s.max.linl = zeros(n1,1);
   s.max.linu(~ieq) = ilam_max.ineqlin(1:n1).*dA(~ieq);
   s.max.linl(~ieq) = ilam_max.ineqlin(n1+1:2*n1).*dA(~ieq);
else
   s.max.linu = [];
   s.max.linl = [];
end
QOISensitivity.Inner = s;
QOISensitivity.Outer = [];
xOpt(:,1) = xmin;
xOpt(:,2) = xmax;

QOIrange.min = Qmin;
QOIrange.max = Qmax;




   function [y,gy] = funxmin(x)
      y = QOIobj.eval(x(id0)');
      gy = zeros(n_variable,1);
      gy(id0) = fg0(x(id0));
   end

   function [y,gy] = funxmax(x)
      y = -QOIobj.eval(x(id0)');
      gy = zeros(n_variable,1);
      gy(id0) = -fg0(x(id0));
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable,2*n_units);
      for i = 1:n_units
         id = idall{i};
         c(2*i-1) = units(i).SurrogateModel.eval(x(id)');
         g(id,2*i-1) = fg{i}(x(id));
      end
      c(2:2:end) = -c(1:2:end);
      g(:,2:2:end) = -g(:,1:2:end);
      c(1:2:end) = c(1:2:end)-UB;
      c(2:2:end) = c(2:2:end)+LB;
      ceq = [];
      geq = [];
   end

end