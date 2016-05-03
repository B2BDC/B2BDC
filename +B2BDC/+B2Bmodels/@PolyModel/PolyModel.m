classdef PolyModel < B2BDC.B2Bmodels.Model
   % A polynomial model object inherited from B2BDC.B2Bmodel.Model
   % superclass
   
   % Created: Nov 11, 2015     Wenyu Li
   
   
   properties
      Basis = []; % A PolyBasis object
      Coefficient = []; % The coefficient for each monomial in the model
   end
   
   properties (SetAccess = private)
      PredictionRange = [];
   end
   
   methods
      function  obj = PolyModel(basis, coefVec, vars, dataStats, err)
         % Create a polynomial model object.
         % The input arguments are:
         %   basis - The PolyBasis object specifies the monomials
         %   coefVec - The coefficient vector for the monomials
         %   vars - A VariableList object contains information of all variables associated with the model 
         %   dataStats - Statistical information about the data (if any) used to fit the model
         %   err - Statistical information about the fitting error (if any)
         if nargin > 0
            if isa(basis,'B2BDC.B2Bvariables.PolyBasis')
               obj.Basis = basis;
            else
               error('Wrong input polynomial basis object')
            end
         else
            obj.Basis = B2BDC.B2Bvariables.PolyBasis();
         end
         if nargin > 1
            if isvector(coefVec) && length(coefVec) == basis.Length
               if isrow(coefVec)
                  obj.Coefficient = coefVec';
               else
                  obj.Coefficient = coefVec;
               end
               id = find(obj.Coefficient == 0);
               if ~isempty(id)
                  obj.Basis.Value(id,:) = [];
                  obj.Coefficient(id) = [];
               end
            else
               error('Wrong input coefficient size or class')
            end
         end
         if nargin > 2
            if vars.Length == basis.Dimension
               obj.Variables = vars;
            else
               error('Input VariablesList object has a wrong dimension')
            end
         end
         if nargin > 3
            obj.DataStats.mx = dataStats.mx;
            obj.DataStats.sx = dataStats.sx;
            obj.DataStats.my = dataStats.my;
            obj.DataStats.sy = dataStats.sy;
         end
         if nargin > 4
            obj.ErrorStats.absMax = err.absMax;
            obj.ErrorStats.absAvg = err.absAvg;
            obj.ErrorStats.relMax = err.relMax;
            obj.ErrorStats.relAvg = err.relAvg;
         end           
      end
      
      function y = eval(obj, X, varObj)
         % Evaluate values of the model at sample points from X. X   
         % should be a tall matrix with size nSample-by-nVariable. 
         % It returns a nSample-by-1 vector of calculated values. If there
         % are 3 input arguments, then the columns of X should be specified
         % by the input VariableList.
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if obj.Basis.Dimension ~= size(X,2)
            error('Wrong input dimension of variables')
         else
            y = zeros(size(X,1),1);
            if ~isempty(obj.DataStats.mx)
               ns = size(X,1);
               mx = obj.DataStats.mx;
               sx = obj.DataStats.sx;
               X = (X-repmat(mx,ns,1))*diag(1./sx);
            end
            for i = 1:size(X,1)
               x = X(i,:);
               for j = 1:obj.Basis.Length
                  bs = obj.Basis.Value(j,:);
                  y(i) = y(i)+obj.Coefficient(j)*prod(x.^bs);
               end
            end
            if ~isempty(obj.DataStats.my);
               my = obj.DataStats.my;
               sy = obj.DataStats.sy;
               y = y*sy+my;
            end
         end
      end
      
      function newPModel = plus(obj,pModel)
         % Returns a new polynomial model which is the sum of the two input
         % component polynomial models.
         if isa(pModel,'B2BDC.B2Bmodels.PolyModel')
            basis1 = obj.Basis;
            basis2 = pModel.Basis;
            var1 = obj.Variables;
            var2 = pModel.Variables;
            newVar = var1.addList(var2);
            basis1 = expandDimension(basis1,var1,newVar);
            basis2 = expandDimension(basis2,var2,newVar);
            newBasis = basis1+basis2;
            coef = zeros(newBasis.Length,1);
            coef1 = obj.Coefficient;
            coef2 = pModel.Coefficient;
            [~,id1,id2] = intersect(basis1.Value,newBasis.Value,'rows');
            coef(id2) = coef(id2)+coef1(id1);
            [~,id1,id2] = intersect(basis2.Value,newBasis.Value,'rows');
            coef(id2) = coef(id2)+coef2(id1);
            newPModel = B2BDC.B2Bmodels.PolyModel(newBasis,coef,newVar);
         else
            error('Wrong input class type')
         end
      end
      
      function newPModel = mtimes(obj, pModel)
         % Returns a new polynomial model which is the product of the two 
         % input component polynomial models.
         if isa(pModel,'B2BDC.B2Bmodels.PolyModel')
            basis1 = obj.Basis;
            basis2 = pModel.Basis;
            coef1 = obj.Coefficient;
            coef2 = pModel.Coefficient;
            var1 = obj.Variables;
            var2 = pModel.Variables;
            newVar = var1.addList(var2);
            basis1 = expandDimension(basis1,var1,newVar);
            basis2 = expandDimension(basis2,var2,newVar);
            newBasis = basis1*basis2;
            coef = zeros(newBasis.Length,1);
            for i = 1:basis1.Length
               for j = 1:basis2.Length
                  tepBasis = basis1.Value(i,:)+ basis2.Value(j,:);
                  tepCoef = coef1(i)*coef2(j);
                  [~,id] = intersect(newBasis.Value,tepBasis,'rows');
                  coef(id) = coef(id)+tepCoef;
               end
            end
            newPModel = B2BDC.B2Bmodels.PolyModel(newBasis,coef,newVar);
         else
            error('Wrong input class type')
         end
      end

      function polyModel = clone(obj)
         % Returns a cloned input polynomial model object.
         basis = obj.Basis;
         Coef = obj.Coefficient;
         vars = obj.Variables;
         dataStat = obj.DataStats;
         err = obj.ErrorStats;
         polyModel = B2BDC.B2Bmodels.PolyModel(basis,Coef,vars,dataStat,err);
      end
      
   end
   
   methods (Hidden = true)
      function dispPoly(obj)
         % Displys the information of the polynomial model.
         n1 = obj.Basis.Length;
         n2 = obj.Basis.Dimension;
         disp(['The polynomial has ' num2str(n2) ' variables and ' num2str(n1) ' monomials'])
         disp('The variables are ordered as')
         disp(obj.VarNames)
         disp('The monomials are ordered as')
         disp(obj.Basis.Value)
         disp('The coefficient are ordered as')
         disp(obj.Coefficient)
      end
   end
   
end

