function [QOIrange, xOpt, QOISensitivity] = predictQOI(obj, QOIobj, B2Bopt)
% Returns the inner and outer bounds, the feasible points of maximum and minimum 
% values of the QOI model subject to the constraints that both dataset units and
% variables are within their respective bounds. The sensitivity of the
% outer bounds with respect to experimental bounds and variable bounds is
% also returned.
% The input arguments are:
%    QOIobj - A B2BDC.B2Bmodel.Model object that defines the QOI
%    B2Bopt - B2BDC.Option object (optional)
% The output are:
%    QOIrange - A structure defines the range of QOI subject to the dataset
%    Xopt - A nVar-by-2 double matrix 

% Created: June 15, 2015   Wenyu Li
%  Modified: July 5, 2015   Wenyu Li (Sensitivity added)
%  Modified: July 14, 2015   Wenyu Li  (Gradient and hessian added)

xOpt = zeros(obj.Variables.Length,2);
if nargin < 3
   B2Bopt = B2BDC.Option();
elseif ~isa(B2Bopt,'B2BDC.Option')
   error('Wrong option object')
end
if ~obj.isConsistent(B2Bopt)
   error('The dataset is inconsistent')
end
vars = obj.Variables.Values;
% LB = [vars.LowerBound]';
% UB = [vars.UpperBound]';
n_variable = length(vars);
LB = -ones(n_variable,1);
UB = -LB;
dx = UB-LB;
units = obj.DatasetUnits.Values;
abE = zeros(length(units),1);
if B2Bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end
flag = quadratictest(units);
if flag && isa(QOIobj,'B2BDC.B2Bmodels.QModel')
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','Hessian',...
      'user-supplied','HessFcn',@hessianfcn,'TolFun',1e-10,'TolCon',1e-10);
else
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','TolFun',1e-10,'TolCon',1e-10);
end
model0 = QOIobj;
n_units = length(units);
allVarnames = obj.VarNames;
[~,id01,id03] = intersect(allVarnames,model0.VarNames);
id02 = [1;id03+1];
if isa(model0,'B2BDC.B2Bmodels.RQModel')
   N0 = zeros(length(id02),length(id02));
   D0 = zeros(length(id02),length(id02));
   N1 = model0.Numerator;
   D1 = model0.Denominator;
   for i0 = 1:length(id02)
      for j0 = 1:length(id02)
         N0(i0,j0) = N1(id02(i0),id02(j0));
         D0(i0,j0) = D1(id02(i0),id02(j0));
      end
   end
elseif isa(model0,'B2BDC.B2Bmodels.QModel')
   Coef0 = zeros(length(id02),length(id02));
   Coef1 = model0.CoefMatrix;
   for i0 = 1:length(id02)
      for j0 = 1:length(id02)
         Coef0(i0,j0) = Coef1(id02(i0),id02(j0));
      end
   end
end

if B2Bopt.Display
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
n_inner = 1;
if n_inner > 1
   minin_result = zeros(n_inner,1);
   xmin = cell(n_inner,1);
   maxin_result = zeros(n_inner,1);
   xmax = cell(n_inner,1);
   if B2Bopt.Display
      h = waitbar(0,'Start to calculate');
   end
   for j = 1:n_inner
      [xmin{j},minin_result(j)] = fmincon(@funxmin,LB+lhsdesign(n_variable,1).*dx,[],[],[],...
         [],LB,UB,@neq,opt);
      [xmax{j},maxin_result(j)] = fmincon(@funxmax,LB+lhsdesign(n_variable,1).*dx,[],[],[],...
         [],LB,UB,@neq,opt);
      if B2Bopt.Display
         waitbar(j/n_inner,h,['Progress: ' int2str(j) ' / ' int2str(n_inner)])
      end
      maxin_result(j) = -maxin_result(j);
   end
   if B2Bopt.Display
      delete(h)
   end
   mx = 0.5*([vars.LowerBound]'+[vars.UpperBound]');
   dx = 0.5*([vars.UpperBound]'-[vars.LowerBound]');
   T = [1, zeros(1,n_variable);
      mx, diag(dx)];
   [minin,id] = min(minin_result);
   xtmp = T*[1;xmin{id}];
   xOpt(:,1) = xtmp(2:end);
   [maxin,id] = max(maxin_result);
   xtmp = T*[1;xmax{id}];
   xOpt(:,2) = xtmp(2:end);
else
   [xmin,minin_result] = fmincon(@funxmin,LB+lhsdesign(n_variable,1).*dx,[],[],[],...
      [],LB,UB,@neq,opt);
   [xmax,maxin_result] = fmincon(@funxmax,LB+lhsdesign(n_variable,1).*dx,[],[],[],...
      [],LB,UB,@neq,opt);
   maxin_result = -maxin_result;
   minin = minin_result;
   maxin = maxin_result;
   mx = 0.5*([vars.LowerBound]'+[vars.UpperBound]');
   dx = 0.5*([vars.UpperBound]'-[vars.LowerBound]');
   T = [1, zeros(1,n_variable);
      mx, diag(dx)];
   xtmp = T*[1;xmin];
   xOpt(:,1) = xtmp(2:end);
   xtmp = T*[1;xmax];
   xOpt(:,2) = xtmp(2:end);
end


if B2Bopt.Display
   disp('=======================================================');
   disp('Calculating outer bound...');
   disp('=======================================================');
end
frac = B2Bopt.ExtraLinFraction;
[minout,minSens] = obj.sedumiminouterbound(QOIobj,frac,abE);
[maxout,maxSens] = obj.sedumimaxouterbound(QOIobj,frac,abE);
%[minout,minSens] = obj.cvxminouterbound(QOIobj,frac);
%[maxout,maxSens] = obj.cvxmaxouterbound(QOIobj,frac);
my = QOIobj.yScale.my;
dy = QOIobj.yScale.dy;
QOIrange.min = dy*[minout, minin]+my;
QOIrange.max = dy*[maxin, maxout]+my;
minSens.expu = dy*minSens.expu;
minSens.expl = dy*minSens.expl;
minSens.varu = dy*minSens.varu;
minSens.varl = dy*minSens.varl;
maxSens.expu = dy*maxSens.expu;
maxSens.expl = dy*maxSens.expl;
maxSens.varu = dy*maxSens.varu;
maxSens.varl = dy*maxSens.varl;
QOISensitivity.min = minSens;
QOISensitivity.max = maxSens;
if B2Bopt.Display
   disp('The calculation is done')
   disp(['Minimum value of QOI is within: [' num2str(minout) ' ' num2str(minin) ']'])
   disp(['Maximum value of QOI is within: [' num2str(maxin) ' ' num2str(maxout) ']'])
end





   function [y,gy] = funxmin(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = ([1;x(id01)]'*N0*[1;x(id01)])/([1;x(id01)]'*D0*[1;x(id01)]);
         gy(id01) = 2*((N0(2:end,2:end)*x(id01)+N0(2:end,1))*([1;x(id01)]'*D0*[1;x(id01)]) - ...
               (D0(2:end,2:end)*x(id01)+D0(2:end,1))*([1;x(id01)]'*N0*[1;x(id01)]))/...
               ([1;x(id01)]'*D0*[1;x(id01)])/([1;x(id01)]'*D0*[1;x(id01)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = [1;x(id01)]'*Coef0*[1;x(id01)];
         gy(id01) = 2*Coef0(2:end,2:end)*x(id01)+2*Coef0(2:end,1);
      end
   end

   function [y,gy] = funxmax(x)
      gy = zeros(n_variable,1);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         y = -([1;x(id01)]'*N0*[1;x(id01)])/([1;x(id01)]'*D0*[1;x(id01)]);
         gy(id01) = -2*((N0(2:end,2:end)*x(id01)+N0(2:end,1))*([1;x(id01)]'*D0*[1;x(id01)]) - ...
               (D0(2:end,2:end)*x(id01)+D0(2:end,1))*([1;x(id01)]'*N0*[1;x(id01)]))/...
               ([1;x(id01)]'*D0*[1;x(id01)])/([1;x(id01)]'*D0*[1;x(id01)]);
      elseif isa(model0,'B2BDC.B2Bmodels.QModel')
         y = -[1;x(id01)]'*Coef0*[1;x(id01)];
         gy(id01) = -2*Coef0(2:end,2:end)*x(id01)-2*Coef0(2:end,1);
      end
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable,2*n_units);
      for i = 1:n_units
         my = units(i).SurrogateModel.yScale.my;
         dy = units(i).SurrogateModel.yScale.dy;
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         l = (l-my)/dy;
         u = (u-my)/dy;
         model = units(i).SurrogateModel;
         [~,id1,id2] = intersect(allVarnames,model.VarNames);
         id3 = [1;id2+1];
         if isa(model,'B2BDC.B2Bmodels.RQModel')
            N = zeros(length(id3),length(id3));
            D = zeros(length(id3),length(id3));
            N1 = model.Numerator;
            D1 = model.Denominator;
            for i1 = 1:length(id3)
               for j1 = 1:length(id3)
                  N(i1,j1) = N1(id3(i1),id3(j1));
                  D(i1,j1) = D1(id3(i1),id3(j1));
               end
            end
            c(2*i-1,1) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)])-u;
            c(2*i,1) =  l-([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id1) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            grad2(id1) = -2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
               (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
               ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(model,'B2BDC.B2Bmodels.QModel')
            quadCoef = zeros(length(id3),length(id3));
            quadCoef2 = model.CoefMatrix;
            for i1 = 1:length(id3)
               for j1 = 1:length(id3)
                  quadCoef(i1,j1) = quadCoef2(id3(i1),id3(j1));
               end
            end
            c(2*i-1,1) = [1;x(id1)]'*quadCoef*[1;x(id1)]-u;
            c(2*i,1) =  l-[1;x(id1)]'*quadCoef*[1;x(id1)];
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            grad2(id1) = -2*quadCoef(2:end,2:end)*x(id1)-2*quadCoef(2:end,1);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         end
      end
      ceq = [];
      geq = [];
   end

   function hessian = hessianfcn(x, lambda)
      hessian = zeros(n_variable,n_variable);
      h0 = model0.Hessian;
      for i2 = 1:length(id03)
         for j2 = 1:length(id03)
            hessian(id01(i2),id01(j2)) = h0(id03(i2),id03(j2));
         end
      end
      for i = 1:n_units
         model = units(i).SurrogateModel;
         modelVar = model.VarNames;
         [~,id1,id2] = intersect(allVarnames,modelVar);
         quadHessian = model.Hessian;
         hessianMatrix1 = zeros(n_variable,n_variable);
         hessianMatrix2 = zeros(n_variable,n_variable);
         for i2 = 1:length(id2)
            for j2 = 1:length(id2)
               hessianMatrix1(id1(i2),id1(j2)) = quadHessian(id2(i2),id2(j2));
               hessianMatrix2(id1(i2),id1(j2)) = -quadHessian(id2(i2),id2(j2));
            end
         end
         hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
            lambda.ineqnonlin(2*i)*hessianMatrix2;
      end
   end

end

