classdef PolyModel < B2BDC.B2Bmodels.Model
   % A polynomial model object inherited from B2BDC.B2Bmodel.Model
   % superclass
   
   % Created: Nov 11, 2015     Wenyu Li
   
   
   properties
      SupportMatrix = []; % A support matrix of size nMonomial-by-nVar
      Coefficient = []; % The coefficient for each monomial in the model
   end
   
   properties (Dependent)
      Degree = [];
   end
   
   methods
      function  obj = PolyModel(sMatrix, coefVec, vars, err)
         % Create a polynomial model object.
         % The input arguments are:
         %   sMatrix - The support matrix specifies the monomials
         %   coefVec - The coefficient vector for the monomials
         %   vars - A VariableList object contains information of all variables associated with the model 
         %   err - Statistical information about the fitting error (if any)
         if nargin > 0
            obj.SupportMatrix = sMatrix;
         end
         if nargin > 1
            if isvector(coefVec) && length(coefVec) == size(sMatrix,1)
               if isrow(coefVec)
                  obj.Coefficient = coefVec';
               else
                  obj.Coefficient = coefVec;
               end
               id = find(obj.Coefficient == 0);
               if ~isempty(id)
                  obj.SupportMatrix(id,:) = [];
                  obj.Coefficient(id) = [];
               end
            else
               error('Wrong input coefficient size or class')
            end
         else
            obj.Coefficient = ones(size(sMatrix,1),1);
         end
         if nargin > 2
            if vars.Length == size(sMatrix,2)
               obj.Variables = vars;
            else
               error('Input VariablesList object has a wrong dimension')
            end
         else
            nV = size(sMatrix,2);
            vars = generateVar([],repmat([-1,1],nV,1));
            obj.Variables = vars;
         end
         if nargin > 3
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
         if obj.Variables.Length ~= size(X,2)
            error('Wrong input dimension of variables')
         else
            c1 = obj.Coefficient;
            s1 = obj.SupportMatrix;
            nM = size(s1,1);
            nS = size(X,1);
            y = zeros(nS,1);
            for i = 1:nS
               x = repmat(X(i,:),nM,1);
               b1 = prod(realpow(x,s1),2);
               y(i) = c1'*b1;
            end
         end
      end
      
      function y = get.Degree(obj)
         s1 = obj.SupportMatrix;
         dm = sum(s1,2);
         y = max(dm);
      end
   end
   
   methods (Hidden = true)
      function dispPoly(obj)
         % Displys the information of the polynomial model.
         s1 = obj.SupportMatrix;
         [n1,n2] = size(s1);
         disp(['The polynomial has ' num2str(n2) ' variables and ' num2str(n1) ' monomials'])
         disp('The variables are ordered as')
         disp(obj.VarNames)
         disp('The monomials are ordered as')
         disp(s1)
         disp('The coefficient are ordered as')
         disp(obj.Coefficient)
      end
   end
   
end

