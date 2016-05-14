classdef RQModel < B2BDC.B2Bmodels.Model
  % A RQModel is a rational quadratic model object inherited from
  % the B2BDC.B2Bmodel.Model superclass
   
   % Created: June 17, 2015     Wenyu Li
   % Modified: Nov 1, 2015     Wenyu Li (error stats added)
   
   properties
      NormalizedNumerator = [];  % The normalized symmetric coefficient matrix of the quadratic numerator of size (nVar+1)-by-(nVar+1)
      NormalizedDenominator = []; % The normalized symmetric coefficient matrix of the quadratic denominator of size (nVar+1)-by-(nVar+1)
      K = []; % Conditional parameter defines the upper bound of denominator over variable domain
   end
   
   properties (Dependent)
      Numerator = []; % The scaled symmetric coefficient matrix of the quadratic numerator of size (nVar+1)-by-(nVar+1)
      Denominator = []; % The scaled symmetric coefficient matrix of the quadratic denominator of size (nVar+1)-by-(nVar+1)
   end
   
   properties (Hidden = true)
      yScale = [];
   end
     
   methods
      function obj = RQModel(Npoly,Dpoly,vars,k,yScale,err)
         % Constructor of the B2BDC.B2Bmodels.RQModel object.
         %
         % The input arguments are:
         %   Npoly - The coefficient matrix of the quadratic numerator
         %   Dpoly - The coefficient matrix of the quadratic denominator
         %   vars - VariableList object containing all model variables
         %   k - The condition parameter defining the upper bound of the denominator over variable ranges
         %   dataStats - Mean and standard deviation of the data used to fit the RQModel
         %   err - Mean and standard deviation of the fitting error of the RQModel
         
         if nargin > 0
            if ~issymmetric(Npoly) || ~issymmetric(Dpoly)
               error('Numerator matrix and denominator matrix must be symmetric')
            else
               obj.NormalizedNumerator = Npoly;
               obj.NormalizedDenominator = Dpoly;
            end 
         end
         if nargin > 2
            obj.Variables = vars;
         end
         if nargin > 3
            if isscalar(k) && k > 1
               obj.K = k;
            else
               error('Parameter K must be a scalar greater than 1')
            end
         end
         if nargin > 4
            if isempty(yScale)
               nVar = obj.Variables.Length;
               x = obj.Variables.makeLHSsample(2*nVar*(nVar+1));
               y = obj.eval(x);
               my = mean(y);
               dy = 0.5*(max(y)-min(y));
               if abs(dy) < 1e-5
                  dy = 1;
               end
               obj.NormalizedNumerator = obj.NormalizedNumerator - my*obj.NormalizedDenominator;
               obj.NormalizedNumerator = obj.NormalizedNumerator/dy;
            else
               my = yScale.my;
               dy = yScale.dy;
            end
            obj.yScale.my = my;
            obj.yScale.dy = dy;
         else
            nVar = obj.Variables.Length;
            x = obj.Variables.makeLHSsample(2*nVar*(nVar+1));
            y = obj.eval(x);
            my = mean(y);
            dy = 0.5*(max(y)-min(y));
            if abs(dy) < 1e-5
               dy = 1;
            end
            obj.yScale.my = my;
            obj.yScale.dy = dy;
            obj.NormalizedNumerator = obj.NormalizedNumerator - my*obj.NormalizedDenominator;
            obj.NormalizedNumerator = obj.NormalizedNumerator/dy;
         end
         if nargin > 5
            obj.ErrorStats.absMax = err.absMax;
            obj.ErrorStats.absAvg = err.absAvg;
            obj.ErrorStats.relMax = err.relMax;
            obj.ErrorStats.relAvg = err.relAvg;
         end  
      end
      
      function y = eval(obj, X, varObj)
         %   Y = EVAL(OBJ, X) evaluates a rational quadratic model at X
         %   sampled points to produce a column vector Y of the
         %   model's output. X is a matrix of size nSample-by-nVariable,
         %   where nSample is the number of samples to be evaluated and
         %   nVariable is the number of variables in the model.
         %
         %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the rational quadratic
         %   model at X sampled points, where X can be of size greater than
         %   nVariable. VAROBJ will specify which columns of X to be used
         %   in evaluating the model to produce a column vector Y.
         
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if size(obj.NormalizedNumerator) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
            x1 = [ones(size(X,1),1), X];
            T = obj.Variables.TransMatrix;
            x2 = T*x1';
            x2 = x2';
            xNew = B2BDC.Fitting.expandBasis(x2(:,2:end));
            Nvec = B2BDC.Fitting.coef2vec(obj.NormalizedNumerator);
            Dvec = B2BDC.Fitting.coef2vec(obj.NormalizedDenominator);
            num = xNew * Nvec;
            den = xNew * Dvec;
            y = num./den;
            if ~isempty(obj.yScale);
               my = obj.yScale.my;
               dy = obj.yScale.dy;
               y = y*dy+my;
            end
         end
      end
      
      function y = get.Numerator(obj)
         T = obj.Variables.TransMatrix;
         Q = obj.NormalizedNumerator;
         y = T'*Q*T;
         y = y*obj.yScale.dy;
         y = y + obj.yScale.my*obj.Denominator;
         y = 0.5*(y + y');
      end
      
      function y = get.Denominator(obj)
         T = obj.Variables.TransMatrix;
         Q = obj.NormalizedDenominator;
         y = T'*Q*T;
      end
         
   end
   
end

