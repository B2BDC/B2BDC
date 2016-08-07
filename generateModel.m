function modelObj = generateModel(varargin)
% Returns a B2BDC.B2Bmodels.Model object. It can take in
% variant input combinations listed below: 
% -------------------------------------------------------------------------
% (a) modelObj = generateModel(GramMatrix, varList)
%  input:
%   GramMatrix - (nVar+1)-by-(nVar+1) symmetric matrix that defines the
%                function of the quadratic model such that
%                M(x) = [1;x]' * GramMatrix * [1;x]
%    varList   - A B2BDC.B2Bvariables.VariableList object contains variable
%                informations (name, domain etc.)
%  output:
%    modelObj  - A B2BDC.B2Bmodels.Qmodel object
%
% (b) modelObj = generateQModel(vecCoef, varList)
%  input:
%    vecCoef   - A (nVar+1)-by-1 vector that defines the function of the
%                linear model such tshat M(x) = vecCoef^T * [1;x]
%    varList   - A B2BDC.B2Bvariables.VariableList object contains variable
%                informations (name, domain etc.)
%  output:
%    modelObj  - A B2BDC.B2Bmodels.Qmodel object
%
% (c) modelObj = generateQModel(GramN, GramD, varList,k)
%  input:
%      GramN   - A (nVar+1)-by-(nVar+1) symmetric matrix that defines the
%                function of the numinator of the rational quadratic model 
%                such that N(x) = [1;x]' * GramN * [1;x]
%      GramD   - A (nVar+1)-by-(nVar+1) symmetric matrix that defines the
%                function of the denominator of the rational quadratic model 
%                such that D(x) = [1;x]' * GramD * [1;x]
%    varList   - A B2BDC.B2Bvariables.VariableList object contains variable
%                informations (name, domain etc.)
%        k     - An optional input that upper bounds the denominator over
%                the variable domain such that 1<=D(x)<=k
%  output:
%    modelObj  - A B2BDC.B2Bmodels.RQmodel object
%
% (d) modelObj = generateQModel(QN, QD, varList,k)
%  input:
%      QN      - A B2BDC.B2Bmodels.QModel object that defines the numerator
%                of the rational quadratic model
%      QD      - A B2BDC.B2Bmodels.QModel object that defines the denominator
%                of the rational quadratic model
%    varList   - A B2BDC.B2Bvariables.VariableList object contains variable
%                informations (name, domain etc.)
%        k     - An optional input that upper bounds the denominator over
%                the variable domain such that 1<=D(x)<=k
%  output:
%    modelObj  - A B2BDC.B2Bmodels.RQmodel object

%  Created: Nov 5, 2015    Wenyu Li

nin = nargin;
%if nin < 2
%   error('Not enought input');
if nin == 1
    modelObj = generateModel([0;1], varargin{:});
elseif nin < 3
   varList = checkVariables(varargin{2});
   nVar = varList.Length;
   if isvector(varargin{1})
      vecCoef = varargin{1};
      if length(vecCoef) ~= nVar+1
         error('The coefficient vector does not match the variable length')
      end
      coefMatrix = zeros(length(vecCoef));
      coefMatrix(1,:) = 0.5*vecCoef;
      coefMatrix = coefMatrix+coefMatrix';
      T = varList.TransMatrix;
      T = inv(T);
      coefMatrix = T'*coefMatrix*T;
      coefMatrix = 0.5*(coefMatrix+coefMatrix');
      if isempty(varList.ExtraQuaConstraint)
         yscale = checkrange;
         my = yscale.my;
         dy = yscale.dy;
         coefMatrix(1,1) = coefMatrix(1,1) - my;
         coefMatrix = coefMatrix / dy;
         modelObj = B2BDC.B2Bmodels.QModel(coefMatrix,varList,yscale);
      else
         modelObj = B2BDC.B2Bmodels.QModel(coefMatrix,varList);
      end
   elseif ismatrix(varargin{1})
      coefMatrix = varargin{1};
      T = varList.TransMatrix;
      T = inv(T);
      coefMatrix = T'*coefMatrix*T;
      coefMatrix = 0.5*(coefMatrix+coefMatrix');
      modelObj = B2BDC.B2Bmodels.QModel(coefMatrix,varList);
   else
      error('Wrong input structure')
   end
elseif nin < 4
   varList = checkVariables(varargin{3});
   if isa(varargin{1},'B2BDC.B2Bmodels.QModel')
      QN = varargin{1};
      QD = varargin{2};
      Ncoef = QN.ScaledCoefMatrix;
      Dcoef = QD.ScaledCoefMatrix;
      T = varList.TransMatrix;
      T = inv(T);
      Ncoef = T'*Ncoef*T;
      Dcoef = T'*Dcoef*T;
      Ncoef = 0.5*(Ncoef+Ncoef');
      Dcoef = 0.5*(Dcoef+Dcoef');
      if Dtest(Dcoef,varList)
         modelObj = B2BDC.B2Bmodels.RQModel(Ncoef,Dcoef,varList);
      else
         error('The denominator function is not always positive over the variable domain')
      end
   elseif ismatrix(varargin{1})
      Ncoef = varargin{1};
      Dcoef = varargin{2};
      T = varList.TransMatrix;
      T = inv(T);
      Ncoef = T'*Ncoef*T;
      Dcoef = T'*Dcoef*T;
      Ncoef = 0.5*(Ncoef+Ncoef');
      Dcoef = 0.5*(Dcoef+Dcoef');
      if Dtest(Dcoef,varList)
         modelObj = B2BDC.B2Bmodels.RQModel(Ncoef,Dcoef,varList);
      else
         error('The denominator function is not always positive over the variable domain')
      end
   else
      error('Wrong input structure')
   end
else
   varList = checkVariables(varargin{3});
   k = varargin{4};
   if isa(varargin{1},'B2BDC.B2Bmodels.QModel')
      QN = varargin{1};
      QD = varargin{2};
      Ncoef = QN.ScaledCoefMatrix;
      Dcoef = QD.ScaledCoefMatrix;
      T = varList.TransMatrix;
      T = inv(T);
      Ncoef = T'*Ncoef*T;
      Dcoef = T'*Dcoef*T;
      Ncoef = 0.5*(Ncoef+Ncoef');
      Dcoef = 0.5*(Dcoef+Dcoef');
      modelObj = B2BDC.B2Bmodels.RQModel(Ncoef,Dcoef,varList,k);
   elseif ismatrix(varargin{1})
      Ncoef = varargin{1};
      Dcoef = varargin{2};
      T = varList.TransMatrix;
      T = inv(T);
      Ncoef = T'*Ncoef*T;
      Dcoef = T'*Dcoef*T;
      Ncoef = 0.5*(Ncoef+Ncoef');
      Dcoef = 0.5*(Dcoef+Dcoef');
      modelObj = B2BDC.B2Bmodels.RQModel(Ncoef,Dcoef,varList,k);
   else
      error('Wrong input structure')
   end
end


   function y = Dtest(Dmatrix,varList)
      Dmodel = B2BDC.B2Bmodels.QModel(Dmatrix,varList);
      testDS = B2BDC.B2Bdataset.Dataset('test');
      tmpMatrix = zeros(size(Dmatrix));
      tmpMatrix(1,1) = 1;
      tmpModel = B2BDC.B2Bmodels.QModel(tmpMatrix, varList);
      tmpDSunit = B2BDC.B2Bdataset.DatasetUnit('test',tmpModel,0,2);
      testDS.addDSunit(tmpDSunit);
      opt = generateOpt('Display',false,'ExtraLinFraction',0);
      minOuter = testDS.sedumiminouterbound(Dmodel, opt.ExtraLinFraction, 0);
      ys = Dmodel.yScale;
      minOuter = minOuter*ys.dy + ys.my;
      if minOuter <= 0
         y = false;
      else
         y = true;
      end
   end

    function b = checkVariables(a)
        if isa(a, 'B2BDC.B2Bvariables.ModelVariable')
            b = B2BDC.B2Bvariables.VariableList();
            b = b.add(a);
        elseif isa(a, 'B2BDC.B2Bvariables.VariableList')
            b = a;
        else
            error('invalid input variables')
        end
    end
 
   function yScale = checkrange
      if ~isempty(varList.ExtraLinConstraint.A)
         A0 = varList.ExtraLinConstraint.A;
         LB = varList.ExtraLinConstraint.LB;
         UB = varList.ExtraLinConstraint.UB;
         [nc,nV] = size(A0);
         H = varList.calBound;
         nS = 2*(nc+nV);
         A = zeros(nS, nS+nV);
         A(:,nV) = [A0; -A0; eye(nV); -eye(nV)];
         A(:,nv+1:end) = eye(nS);
         b = [UB; -LB; H(:,2); -H(:,1)];
         k.f = nV;
         k.l = nS;
         pars.fid = 0;
         c1 = vecCoef(2:end);
         tmpx = sedumi(A,b,c1,k,pars);
         miny = c1'*tmpx(1:nV)+c1(1);
         tmpx = sedumi(A,b,-c1,k,pars);
         maxy = -c1'*tmpx(1:nV)+c1(1);
      else
         c = vecCoef;
         miny = c(1);
         maxy = c(1);
         nV = varList.Length;
         H = varList.calBound;
         for i = 1:nV
            if c(i+1) >= 0
               maxy = maxy + H(i,2)*c(i+1);
               miny = miny + H(i,1)*c(i+1);
            else
               maxy = maxy + H(i,1)*c(i+1);
               miny = miny + H(i,2)*c(i+1);
            end
         end
      end
      yScale.my = 0.5*(miny+maxy);
      yScale.dy = 0.5*(maxy-miny);
   end
end


