function [Qunits, Qx, Qextra, nextra, extraIdx, L, idRQ] = getInequalQuad(obj,bds,frac)
% Return the quadratic form of inequality constraints of the
% B2BDC.B2Bdataset.Dataset object. The returned quadratic form matrix is
% with respect to all the active variables of the dataset.
% Input:
%  bds  -  A nUnits-by-2 matrix defines the lower and upper bounds for each
%          surrogate model in the dataset unit.
%  frac -  Fraction of extra constraints used in the optimization
%          0 < frac < 100, if frac == -1, then automatically linear paris
%          with influence factor greater than 5% of the most influential pair
%          will be included
% Output:
%  Qunits - Inequality quadratic form of the observation QOI to be within
%           the uncertainty bounds
%    Qx   - Inequality quadratic form of the linear box constraints of all
%           the active variables of the dataset
%  Qextra - Inequality quadratic form of extra combination pairs of linear
%           constraints on the active variables of the dataset
%  n_extra - Number of extra variable pairs used in the SDP
% extraIdx - A n_extra-by-2 matrix stores the index of each variable pair
% All the inequality quadratic matrix is formed as Q <= 0

% Created: August 10, 2015     Wenyu Li

% eps = 1e-7;
vars = obj.Variables;
A0 = vars.ExtraLinConstraint.A;
if ~isempty(A0)
   lLB = vars.ExtraLinConstraint.LB;
   lUB = vars.ExtraLinConstraint.UB;
end
nlX = size(A0,1);
varName = obj.VarNames;
n_variable = vars.Length;
xbds = vars.calBound;
LB = xbds(:,1);
UB = xbds(:,2);
units = obj.DatasetUnits.Values;
n_units = length(units);
Qunits = cell(2*n_units,1);
Qx = cell(n_variable+nlX,1);
idRQ = false(n_units,1);
for i = 1:n_units
   model = units(i).SurrogateModel;
   modelVar = model.VarNames;
   [~,id1,id2] = intersect(varName,modelVar);
   id1 = [1;id1+1];
   id2 = [1;id2+1];
   if isa(model,'B2BDC.B2Bmodels.RQModel')
      idRQ(i) = true;
      N = model.Numerator;
      D = model.Denominator;
      constrN = zeros(n_variable+1,n_variable+1);
      constrD = zeros(n_variable+1,n_variable+1);
      for j = 1:length(id2)
         for k = 1:length(id2)
            constrN(id1(j),id1(k)) = N(id2(j),id2(k));
            constrD(id1(j),id1(k)) = D(id2(j),id2(k));
         end
      end
      Qunits{2*i-1} = constrN-bds(i,2)*constrD;
      Qunits{2*i} = bds(i,1)*constrD-constrN;
   elseif isa(model,'B2BDC.B2Bmodels.QModel')
      Coef = model.CoefMatrix;
      constrMatrix1 = zeros(n_variable+1,n_variable+1);
      constrMatrix2 = zeros(n_variable+1,n_variable+1);
      for j = 1:length(id2)
         for k = 1:length(id2)
            constrMatrix1(id1(j),id1(k)) = Coef(id2(j),id2(k));
            constrMatrix2(id1(j),id1(k)) = -Coef(id2(j),id2(k));
         end
      end
      constrMatrix1(1,1) = constrMatrix1(1,1)-bds(i,2);
      constrMatrix2(1,1) = constrMatrix2(1,1)+bds(i,1);
      Qunits{2*i-1} = constrMatrix1;
      Qunits{2*i} = constrMatrix2;
   end
end
for i = 1:n_variable
   lx = LB(i);
   ux = UB(i);
   Qx{i} = sparse([1,1,i+1,i+1],[1,i+1,1,i+1],...
      [lx*ux,-(lx+ux)/2,-(lx+ux)/2,1],n_variable+1,n_variable+1);
end
for i = 1:nlX
   ai = A0(i,:);
   E = zeros(n_variable+1);
   E(2:end,2:end) = ai'*ai;
   E(1,1) = lLB(i)*lUB(i);
   E(1,2:end) = -0.5*(lLB(i)+lUB(i))*ai;
   E(2:end,1) = -0.5*(lLB(i)+lUB(i))*ai';
   Qx{n_variable+i} = E;
end
if ~isempty(A0)
   Le = A0;
end
Lv = speye(n_variable);
if ~isempty(A0)
   L = [Lv;Le];
else
   L = Lv;
end
if frac == -1
   [Qextra, extraIdx] = getExtraQ;
%    tic;
%    [Qextra, extraIdx] = getExtraQ;
%    [~, Idx] = B2BDC.B2Bdataset.Dataset.findConicHull([Qunits;Qx],Qe);
%    Qextra = Qe(Idx);
%    toc;
%    extraIdx = extraIdx(Idx,:);
   nextra = length(Qextra);
elseif frac == 0
   Qextra = {};
   extraIdx = [];
   nextra = 0;
elseif frac == -2
   load idx
   [Qe, extraIdx] = getExtraQ;
   Qextra = Qe(Idx);
   extraIdx = extraIdx(Idx,:);
   nextra = length(Qextra);
else
   [Qe, extraIdx] = getExtraQ;
   nextra = round(length(Qe)*frac*0.01);
   [~, Idx] = B2BDC.B2Bdataset.Dataset.approxConicHull([Qunits;Qx],Qe,nextra);
   Qextra = Qe(Idx);
   extraIdx = extraIdx(Idx,:);
end




   function [Qe, idx] = getExtraQ
      if ~isempty(A0)
         lb = [LB; lLB];
         ub = [UB; lUB];
      else
         lb = LB;
         ub = UB;
      end
      [nL,nVar] = size(L);
      n1 = nchoosek(nL,2)*4;
      Qe = cell(n1,1);
      idx = zeros(n1,3);
      count = 1;
      for i1 = 1:nL-1
         for j1 = i1+1:nL
            L1 = L(i1,:);
            L2 = L(j1,:);
            l1 = lb(i1);
            l2 = lb(j1);
            u1 = ub(i1);
            u2 = ub(j1);
            %lb lb
            Q = zeros(nVar+1);
            Q(1,1) = -0.5*l1*l2;
            Q(1,2:end) = 0.5*(l1*L2+l2*L1);
            Q(2:end,2:end) = -0.5*L1' * L2;
            Q = Q+Q';
            Qe{count} = Q;
            idx(count,:) = [i1,j1,1];
            count = count+1;
            %lb ub
            Q = zeros(nVar+1);
            Q(1,1) = 0.5*l1*u2;
            Q(1,2:end) = -0.5*(l1*L2+u2*L1);
            Q(2:end,2:end) = 0.5*L1' * L2;
            Q = Q+Q';
            Qe{count} = Q;
            idx(count,:) = [i1,j1,2];
            count = count+1;
            %ub lb
            Q = zeros(nVar+1);
            Q(1,1) = 0.5*u1*l2;
            Q(1,2:end) = -0.5*(u1*L2+l2*L1);
            Q(2:end,2:end) = 0.5*L1' * L2;
            Q = Q+Q';
            Qe{count} = Q;
            idx(count,:) = [i1,j1,3];
            count = count+1;
            %ub ub
            Q = zeros(nVar+1);
            Q(1,1) = -0.5*u1*u2;
            Q(1,2:end) = 0.5*(u1*L2+u2*L1);
            Q(2:end,2:end) = -0.5*L1' * L2;
            Q = Q+Q';
            Qe{count} = Q;
            idx(count,:) = [i1,j1,4];
            count = count+1;
         end
      end
   end

end