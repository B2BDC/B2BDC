classdef QModel < B2BDC.B2Bmodels.Model
   % A QModel is a quadratic model object inherited from
   % B2BDC.B2Bmodel.Model superclass
   
   % Created: June 17, 2015     Wenyu Li
   % Modified: Nov 1, 2015     Wenyu Li (error stats added)
   
   properties (Dependent)
      NormalizedCoefMatrix = []; % The normalized symmetric coefficient matrix of the quadratic model of size (nVar+1)-by-(nVar+1)
   end
   
   properties 
      CoefMatrix = [];  % The original symmetric coefficient matrix of the quadratic model of size (nVar+1)-by-(nVar+1)
   end
   
   properties (Dependent, Hidden = true)
      Hessian = [];
   end
   
   properties (Hidden = true)
      yScale = [];
   end
   
   methods
      function obj = QModel(coef,vars,yScale,err)
            % Constructor of the B2BDC.B2Bmodels.QModel object.
            %
            % The input arguments are:
            %   coef - The coefficient matrix of the quadratic model
            %   vars - VariableList object containing all model variables
            %   dataStats - Mean and standard deviation of the data used to fit the QModel
            %   err - Mean and standard deviation of the fitting error of the QModel
            
         if nargin > 0
            if issymmetric(coef)
               obj.CoefMatrix = coef;
            else
               error('Coefficient matrix must be symmetric')
            end
         end
         if nargin > 1
            obj.Variables = vars;
         end
         if nargin > 2
            if isempty(yScale)
               nVar = obj.Variables.Length;
               x = obj.Variables.makeLHSsample(nVar*(nVar+1));
               y = obj.eval(x);
               my = mean(y);
               dy = 0.5*(max(y)-min(y));
               if abs(dy) < 1e-3
                  dy = 1;
               end
               yScale.my = my;
               yScale.dy = dy; 
            end
            obj.yScale = yScale;
         else
            nVar = obj.Variables.Length;
            x = obj.Variables.makeLHSsample(nVar*(nVar+1));
            y = obj.eval(x);
            my = mean(y);
            dy = 0.5*(max(y)-min(y));
            if abs(dy) < 1e-3
               dy = 1;
            end
            obj.yScale.my = my;
            obj.yScale.dy = dy;
         end
         if nargin > 3
            obj.ErrorStats.absMax = err.absMax;
            obj.ErrorStats.absAvg = err.absAvg;
            obj.ErrorStats.relMax = err.relMax;
            obj.ErrorStats.relAvg = err.relAvg;
         end  
      end
      
      function y = eval(obj, X, varObj)
        %   Y = EVAL(OBJ, X) evaluates a quadratic model at X sampled
        %   points to produce a column vector Y of the model's output. X
        %   is a matrix of size nSample-by-nVariable, where nSample is the
        %   number of samples to be evaluated and nVariable is the number
        %   of variables in the model.
        %
        %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the quadratic model at
        %   X sampled points, where X can be of size greater than nVariable.
        %   VAROBJ will specify which columns of X to be used in evaluating
        %   the model to produce a column vector Y of the model's output.
        
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if size(obj.CoefMatrix,1) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
            nSample = size(X,1);
            x1 = [ones(nSample,1), X];
            xNew = B2BDC.Fitting.expandBasis(x1(:,2:end));
            covec = B2BDC.Fitting.coef2vec(obj.CoefMatrix);
            y = xNew*covec;
%             coef = obj.CoefMatrix;
%             y = diag(x1 * coef * x1');
         end
      end
      
      function y = get.Hessian(obj)
         y = 2 * obj.CoefMatrix(2:end,2:end);
      end
      
      function y = get.NormalizedCoefMatrix(obj)
         T = obj.Variables.TransMatrix;
         Q = obj.CoefMatrix;
         y = T'\Q/T;
         y(1,1) = y(1,1) - obj.yScale.my;
         y = y/obj.yScale.dy;
         y = 0.5*(y + y');
      end
   end
   
   methods (Hidden = true)
      function [A,b] = calGradient(obj)
         quadCoef = obj.CoefMatrix(2:end,2:end);
         A = 2 * quadCoef;
         b = 2 * obj.CoefMatrix(2:end,1);
      end
   end
   
end

