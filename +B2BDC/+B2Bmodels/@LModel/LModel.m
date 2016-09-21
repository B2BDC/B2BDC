classdef LModel < B2BDC.B2Bmodels.Model
   % A QModel is a quadratic model object inherited from
   % B2BDC.B2Bmodel.Model superclass
   
   % Created: Sep 20, 2016     Wenyu Li
   
   properties (Dependent)
      NormalizedCoefVec = []; % The normalized coefficient vector of the linear model
   end
   
   properties
      CoefVec = []; % The original coefficient vector of the linear model
   end
   
   properties (Hidden = true)
      yScale = [];
   end
   
   methods
      function obj = LModel(coef,vars,yScale,err)
         % Constructor of the B2BDC.B2Bmodels.LModel object.
            %
            % The input arguments are:
            %   coef - The coefficient vector of the linear model
            %   vars - VariableList object containing all model variables
            %   yScale - Mean and standard deviation of the data used to fit the QModel
            %   err - Mean and standard deviation of the fitting error of the QModel
            
            if nargin > 0
               if size(coef,1) == 1
                  obj.CoefVec = coef';
               else
                  obj.CoefVec = coef;
               end
            end
            if nargin > 1
               if length(obj.CoefVec) ~= vars.Length+1
                  error('Mismatched variable and coefficient dimension');
               else
                  obj.Variables = vars;
               end
            end
            if nargin > 2
               if isempty(yScale)
                  c = coef;
                  miny = c(1);
                  maxy = c(1);
                  nV = vars.Length;
                  H = vars.calBound;
                  for i = 1:nV
                     if c(i+1) >= 0
                        maxy = maxy + H(i,2)*c(i+1);
                        miny = miny + H(i,1)*c(i+1);
                     else
                        maxy = maxy + H(i,1)*c(i+1);
                        miny = miny + H(i,2)*c(i+1);
                     end
                  end
                  yScale.my = 0.5*(miny+maxy);
                  yScale.dy = 0.5*(maxy-miny);
                  obj.yScale = yScale;
               else
                  obj.yScale = yScale;
               end
            else
               c = coef;
               miny = c(1);
               maxy = c(1);
               nV = vars.Length;
               H = vars.calBound;
               for i = 1:nV
                  if c(i+1) >= 0
                     maxy = maxy + H(i,2)*c(i+1);
                     miny = miny + H(i,1)*c(i+1);
                  else
                     maxy = maxy + H(i,1)*c(i+1);
                     miny = miny + H(i,2)*c(i+1);
                  end
               end
               yScale.my = 0.5*(miny+maxy);
               yScale.dy = 0.5*(maxy-miny);
               obj.yScale = yScale;
            end
            if nargin > 3
               obj.ErrorStats.absMax = err.absMax;
               obj.ErrorStats.absAvg = err.absAvg;
               obj.ErrorStats.relMax = err.relMax;
               obj.ErrorStats.relAvg = err.relAvg;
            end     
      end
      
      function y = eval(obj, X, varObj)
         %   Y = EVAL(OBJ, X) evaluates a linear model at X sampled
         %   points to produce a column vector Y of the model's output. X
         %   is a matrix of size nSample-by-nVariable, where nSample is the
         %   number of samples to be evaluated and nVariable is the number
         %   of variables in the model.
         %
         %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the linear model at
         %   X sampled points, where X can be of size greater than nVariable.
         %   VAROBJ will specify which columns of X to be used in evaluating
         %   the model to produce a column vector Y of the model's output.
         
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if length(obj.CoefVec) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
            nSample = size(X,1);
            x1 = [ones(nSample,1), X];
            y = x1*obj.CoefVec;
         end
      end
      
      function y = get.NormalizedCoefVec(obj)
         T = obj.Variables.TransMatrix;
         L = obj.CoefVec;
         y = T'\L;
      end
   end
   
end

