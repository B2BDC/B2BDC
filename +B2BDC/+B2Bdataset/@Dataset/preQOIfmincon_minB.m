function [Qmin, Qmax, ss, xOpt, alpha] = preQOIfmincon_minB(obj,q,disflag,rflag,b2bopt,alpha)
% subfunction to calculate QOI inner bound with fmincon

%  Created: Oct 10, 2016     Wenyu Li

if rflag
   error('The predicted QOI has variables not in the dataset');
end
vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
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
   n1 = size(A1,1);
   n2 = size(Aeq,1);
   A1 = [zeros(n1,1), A1];
   if n2 ~= 0
      Aeq = [zeros(n2,1), Aeq];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
xFea = obj.FeasiblePoint;
if obj.ModelDiscrepancyFlag
   MDvar = obj.ModelDiscrepancy.Variables;
   nMD = MDvar.Length;
   MDbd = MDvar.calBound;
   LB = [LB; MDbd(:,1)];
   UB = [UB; MDbd(:,2)];
   if obj.ParameterDiscrepancyFlag
      PDvar = obj.ParameterDiscrepancy.Variables;
      nPD = PDvar.Length;
      PDbd = PDvar.calBound;
      LB = [LB; PDbd(:,1)];
      UB = [UB; PDbd(:,2)];
   else
      nPD = 0;
   end
else
   nMD = 0;
   if obj.ParameterDiscrepancyFlag
      PDvar = obj.ParameterDiscrepancy.Variables;
      nPD = PDvar.Length;
      PDbd = PDvar.calBound;
      LB = [LB; PDbd(:,1)];
      UB = [UB; PDbd(:,2)];
   else
      nPD = 0;
   end
end
if nMD > 0
   xFea = [xFea; obj.ModelDiscrepancy.FeasiblePoint];
end
if nPD > 0
   xFea = [xFea; obj.ParameterDiscrepancy.FeasiblePoint];
end
dx = UB-LB;
[idall,Qall,Nall,Dall,APD,bPD] = obj.getQ_RQ_expansion(q);
id0 = idall{end};
Q0 = Qall{end};
N0 = Nall{end};
D0 = Dall{end};
if ~isempty(A1)
   A1 = [A1 zeros(size(A1,1),nMD+nPD)]; 
end
if ~isempty(Aeq)
   Aeq = [Aeq zeros(size(Aeq,1),nMD+nPD)];
end
if ~isempty(APD)
   infFlag = bPD == inf;
   APD(infFlag,:) = [];
   bPD(infFlag) = [];
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
bds = obj.calBound;
bds = bds(1:n_units,:);
bds = diff(bds,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
opt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point',...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
nStart = b2bopt.LocalStart;
% optimal B
xStart = xFea';
if nStart > 1
   xSample = obj.collectHRsamples(nStart-1,[],b2bopt);
   xStart = [xStart; xSample.x];
end
xmin = cell(nStart,1);
ymin = zeros(nStart,1);
minflag = zeros(nStart,1);
for j = 1:nStart
   [xmin{j},ymin(j),minflag(j)] = fmincon(@minB,xStart(j,:)',...
      [A1;APD],[B1;bPD],Aeq,Beq,LB,UB,@neq,opt);
end
ymin(minflag<=0) = inf;
[aopt,idB] = min(ymin);
xmin = xmin{idB};
Bopt = xmin(n_variable+1:n_variable+nMD);
ida = find(alpha>aopt,1);
if ~isempty(ida)
   alpha = [aopt, alpha(ida:end)];
else
   alpha = aopt;
end
nTest = length(alpha);
% aa = [zeros(nMD,n_variable) eye(nMD) zeros(nMD,nPD)];
% Aeq = [Aeq;aa];
% Beq = [Beq;Bopt];
Qmin = zeros(nTest,1);
Qmax = zeros(nTest,1);
ss = cell(nTest,1);
xOpt = cell(nTest,1);
for k = 1:nTest
   tol = alpha(k);
   xStart = repmat(xmin',2*nStart,1);
   if nStart > 1
      xStart(3:end,[1:n_variable n_variable+nMD+1:n_variable+nMD+nPD]) = ...
         xStart(3:end,[1:n_variable n_variable+nMD+1:n_variable+nMD+nPD]) + 0.05*randn(2*(nStart-1),n_variable+nPD);
   end
   xmin = cell(nStart,1);
   ymin = zeros(nStart,1);
   minflag = zeros(nStart,1);
   ilam_min = cell(nStart,1);
   xmax = cell(nStart,1);
   ymax = zeros(nStart,1);
   maxflag = zeros(nStart,1);
   ilam_max = cell(nStart,1);
   for j = 1:nStart
      [xmin{j},ymin(j),minflag(j),~,ilam_min{j}] = fmincon(@funxmin,xStart(2*j-1,:)',...
         [A1;APD],[B1;bPD],Aeq,Beq,LB,UB,@neq_B,opt);
      [xmax{j},ymax(j),maxflag(j),~,ilam_max{j}] = fmincon(@funxmax,xStart(2*j,:)',...
         [A1;APD],[B1;bPD],Aeq,Beq,LB,UB,@neq_B,opt);
   end
   ymin(minflag<=0) = inf;
   [Qmin(k),ID] = min(ymin);
   xmin = xmin{ID};
   ilam_min = ilam_min{ID};
   ymax(maxflag<=0) = -inf;
   [Qmax(k),ID] = min(ymax);
   Qmax(k) = -Qmax(k);
   xmax = xmax{ID};
   ilam_max = ilam_max{ID};
   
   s.min.expu = -ilam_min.ineqnonlin(1:2:2*n_units).*bds;
   s.min.expl = -ilam_min.ineqnonlin(2:2:2*n_units).*bds;
   s.min.varu = -ilam_min.upper.*dx;
   s.min.varl = -ilam_min.lower.*dx;
   n3 = size(A1,1);
   if ~isempty(A1)
      s.min.linu = zeros(n1+n2,1);
      s.min.linl = zeros(n1+n2,1);
      s.min.linu(~ieq) = -ilam_min.ineqlin(1:n1).*dA(~ieq);
      s.min.linl(~ieq) = -ilam_min.ineqlin(n1+1:2*n1).*dA(~ieq);
   else
      s.min.linu = [];
      s.min.linl = [];
   end
   if ~isempty(APD)
      n4 = 0.5*size(APD,1);
      dAPD = bPD(1:2:end)+bPD(2:2:end);
      s.min.linu = [s.min.linu ; -ilam_min.ineqlin(2*n3+1:2:2*(n3+n4)).*dAPD];
      s.min.linl = [s.min.linl ; -ilam_min.ineqlin(2*n3+2:2:2*(n3+n4)).*dAPD];
   end
   s.max.expu = ilam_max.ineqnonlin(1:2:2*n_units).*bds;
   s.max.expl = ilam_max.ineqnonlin(2:2:2*n_units).*bds;
   s.max.varu = ilam_max.upper.*dx;
   s.max.varl = ilam_max.lower.*dx;
   if ~isempty(A1)
      s.max.linu = zeros(n1+n2,1);
      s.max.linl = zeros(n1+n2,1);
      s.max.linu(~ieq) = ilam_max.ineqlin(1:n1).*dA(~ieq);
      s.max.linl(~ieq) = ilam_max.ineqlin(n1+1:2*n1).*dA(~ieq);
   else
      s.max.linu = [];
      s.max.linl = [];
   end
   if ~isempty(APD)
      n4 = 0.5*size(APD,1);
      dAPD = bPD(1:2:end)+bPD(2:2:end);
      s.max.linu = [s.max.linu ; ilam_max.ineqlin(2*n3+1:2:2*(n3+n4)).*dAPD];
      s.max.linl = [s.max.linl ; ilam_max.ineqlin(2*n3+2:2:2*(n3+n4)).*dAPD];
   end
   ss{k} = s;
   % xOpt.min.x = xmin(1:n_variable);
   % xOpt.max.x = xmax(1:n_variable);
   % if nMD > 0
   %    xOpt.min.xMD = xmin(n_variable+1:n_variable+nMD);
   %    xOpt.max.xMD = xmax(n_variable+1:n_variable+nMD);
   % end
   % if nPD > 0
   %    xOpt.min.xPD = xmin(n_variable+nMD+1:n_variable+nMD+nPD);
   %    xOpt.max.xPD = xmax(n_variable+nMD+1:n_variable+nMD+nPD);
   % end
   xOpt{k}.min.x = xmin;
   xOpt{k}.min.dimension = [n_variable nMD nPD];
   xOpt{k}.max.x = xmax;
   xOpt{k}.max.dimension = [n_variable nMD nPD];
end


   function [y,gy] = funxmin(x)
      gy = zeros(n_variable+nMD+nPD,1);
      if ~isempty(N0)
         y = ([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
         gy(id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif ~isempty(Q0)
         y = [1;x(id0)]'*Q0*[1;x(id0)];
         gy(id0) = 2*(Q0(2:end,2:end)*x(id0)+Q0(2:end,1));
      end
   end

   function [y,gy] = funxmax(x)
      gy = zeros(n_variable+nMD+nPD,1);
      if ~isempty(N0)
         y = -([1;x(id0)]'*N0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
         gy(id0) = -2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif ~isempty(Q0)
         y = -[1;x(id0)]'*Q0*[1;x(id0)];
         gy(id0) = -2*(Q0(2:end,2:end)*x(id0)+Q0(2:end,1));
      end
   end

   function [y,gy] = minB(x)
      y = zeros(1,n_units);
      gy = zeros(n_variable+nMD+nPD,n_units);
      for i = 1:n_units
         id = idall{i};
         id1 = id(id<=n_variable);
         id2 = id(id>n_variable+nMD);
         nn1 = length(id1);
         nn2 = length(id2);
         nn3 = length(id)-nn1-nn2;
         id3 = id(nn1+1:nn1+nn3);
         Coef = Qall{i};
         M = ([1;x([id1;id2])])'*Coef([1;id1+1;id2+1],[1;id1+1;id2+1])*[1;x([id1;id2])];
         B = 2*Coef(1,id3+1)*x(id3);
         y(i) = B/M;
         gy([id1;id2],i) = -2*sign(y(i))*B*(Coef([id1+1;id2+1],[id1+1;id2+1])*x([id1;id2])+Coef([id1+1;id2+1],1))/M^2;
         gy(id3,i) = 2*sign(y(i))*Coef(id3+1,1)/M;
      end
      y = mean(abs(y));
      gy = mean(gy,2);
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+nPD,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)])-u;
            c(2*i,1) =  l-([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad2(id) = -2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            Coef = Qall{i};
            c(2*i-1,1) = [1;x(id)]'*Coef*[1;x(id)]-u;
            c(2*i,1) =  l-[1;x(id)]'*Coef*[1;x(id)];
            g(id,2*i-1) = 2*Coef(2:end,2:end)*x(id)+2*Coef(2:end,1);
            g(id,2*i) = -g(id,2*i-1);
         end
      end
      ceq = [];
      geq = [];
   end

   function [c,ceq,g,geq] = neq_B(x)
      c = zeros(2*n_units+1,1);
      g = zeros(n_variable+nMD+nPD,2*n_units+1);
      y = zeros(1,n_units);
      gy = zeros(n_variable+nMD+nPD,n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         id = idall{i};
         id1 = id(id<=n_variable);
         id2 = id(id>n_variable+nMD);
         nn1 = length(id1);
         nn2 = length(id2);
         nn3 = length(id)-nn1-nn2;
         id3 = id(nn1+1:nn1+nn3); 
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)])-u;
            c(2*i,1) =  l-([1;x(id)]'*N*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad1 = zeros(n_variable,1);
            grad2 = zeros(n_variable,1);
            grad1(id) = 2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            grad2(id) = -2*((N(2:end,2:end)*x(id)+N(2:end,1))*([1;x(id)]'*D*[1;x(id)]) - ...
               (D(2:end,2:end)*x(id)+D(2:end,1))*([1;x(id)]'*N*[1;x(id)]))/...
               ([1;x(id)]'*D*[1;x(id)])/([1;x(id)]'*D*[1;x(id)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            Coef = Qall{i};
            c(2*i-1,1) = [1;x(id)]'*Coef*[1;x(id)]-u;
            c(2*i,1) =  l-[1;x(id)]'*Coef*[1;x(id)];
            g(id,2*i-1) = 2*Coef(2:end,2:end)*x(id)+2*Coef(2:end,1);
            g(id,2*i) = -g(id,2*i-1);
            M = [1;x([id1;id2])]'*Coef([1;id1+1;id2+1],[1;id1+1;id2+1])*[1;x([id1;id2])];
            B = 2*Coef(1,id3+1)*x(id3);
            y(i) = B/M;
            gy([id1;id2],i) = -2*sign(y(i))*B*(Coef([id1+1;id2+1],[id1+1;id2+1])*x([id1;id2])+Coef([id1+1;id2+1],1))/M^2;
            gy(id3,i) = 2*sign(y(i))*Coef(id3+1,1)/M;
         end
      end
      c(end) = mean(abs(y))-tol;
      g(:,end) = mean(gy,2);
      ceq = [];
      geq = [];
   end

   function hessian = hessianfcn(x, lambda)
      hessian = zeros(n_variable);
      if isa(model0,'B2BDC.B2Bmodels.RQModel')
         hessian(id0,id0) = 2*((N0(2:end,2:end)*x(id0)+N0(2:end,1))*([1;x(id0)]'*D0*[1;x(id0)]) - ...
            (D0(2:end,2:end)*x(id0)+D0(2:end,1))*([1;x(id0)]'*N0*[1;x(id0)]))/...
            ([1;x(id0)]'*D0*[1;x(id0)])/([1;x(id0)]'*D0*[1;x(id0)]);
      elseif isa(model0,'B2BDC.B2Bmodels.RQModel')
         hessian(id0,id0) = model0.Hessian;
      end
      for i = 1:n_units
         id = idall{i};
         quadHessian = 2*Qall{i}(2:end,2:end);
         hessianMatrix1 = zeros(n_variable);
         hessianMatrix2 = zeros(n_variable);
         hessianMatrix1(id,id) = quadHessian;
         hessianMatrix2(id,id) = -quadHessian;
         hessian = hessian+lambda.ineqnonlin(2*i-1)*hessianMatrix1+...
            lambda.ineqnonlin(2*i)*hessianMatrix2;
      end
   end

end