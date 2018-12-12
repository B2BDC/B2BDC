classdef NNModel < B2BDC.B2Bmodels.Model
   % A NNmodel is a neural network model object inherited from
   % B2BDC.B2Bmodel.Model superclass
   
   % Created: Feb 20, 2018     Wenyu Li
   
   properties
      Weights = {};  % The weights used in the network layers
      Biases = {};  % The biases used in the network layers
      Constant = 0; % The constant added to the NNmodel output
   end
   
   properties(SetAccess = private)
      ActivationFcn = {}; % The activation function used in the network layers
      ActivationGradientFcn = {}; % The gradient of the activation function used in the network layers
      NetGradientFcn = []; % The gradient of the network output with respect to input
   end
   
   properties(Dependent = true)
      LayerNumber = [];  % The number of layers of the network
   end
   
   properties(Hidden = true, SetAccess = private)
      NNobject = []; % The MATLAB net object
      xTransform = [];  % Pre-processing parameter for input
      yTransform = [];  % Post-processing parameter for input
   end
   
   methods
      function obj = NNModel(nnobj,vars,err)
            % Constructor of the B2BDC.B2Bmodels.NNModel object.
            %
            % The input arguments are:
            %   nnobj - A MATLAB network object
            %   vars - VariableList object containing all model variables
            %   err - Mean and maximum fitting error of the NNModel
            
            if nargin > 0 && isa(nnobj,'network')
               obj.NNobject = nnobj;
               nLayer = nnobj.numLayers;
               weights = cell(nLayer,1);
               weights{1} = nnobj.IW{1};
               lw = nnobj.LW(2:end,1:end-1);
               for i = 2:nLayer
                  weights{i} = lw{i-1,i-1};
               end
               obj.Weights = weights;
               obj.Biases = nnobj.b;
               actfcn = cell(nLayer,1);
               gfcn = cell(nLayer,1);
               for i = 1:nLayer
                  switch nnobj.layers{i}.transferFcn
                     case 'tansig'
                        actfcn{i} = @tansig;
                        gfcn{i} = @(x) 1-tansig(x).^2;
                     case 'logsig'
                        actfcn{i} = @logsig;
                        gfcn{i} = @(x) logsig(x).*(1-logsig(x));
                     case 'purelin'
                        actfcn{i} = @purelin;
                        gfcn{i} = @(x) ones(length(x),1);
                  end
               end
               obj.ActivationFcn = actfcn;
               obj.ActivationGradientFcn = gfcn;
               inp = nnobj.inputs{1}.processSettings{1};
               x.Range = inp.xrange';
               x.Translation = inp.xmin';
               obj.xTransform = x;
               outp = nnobj.outputs{end}.processSettings{1};
               y.Range = outp.xrange';
               y.Translation = outp.xmin';
               obj.yTransform = y;
               obj.NetGradientFcn = calculateGradient(obj);
            else
               error('The input is not a network object');
            end   
            if nargin > 1
               obj.Variables = vars;
               nVar = vars.Length;
               if size(weights{1},2) ~= nVar
                  error('The dimension of the network inputs is not consistent')
               end
            end
            if nargin > 2
               obj.ErrorStats.absMax = err.absMax;
               obj.ErrorStats.absAvg = err.absAvg;
               obj.ErrorStats.relMax = err.relMax;
               obj.ErrorStats.relAvg = err.relAvg;
            end
      end
      
      function y = eval(obj, X, varObj)
        %   Y = EVAL(OBJ, X) evaluates a network model at X sampled
        %   points to produce a column vector Y of the model's output. X
        %   is a matrix of size nSample-by-nVariable, where nSample is the
        %   number of samples to be evaluated and nVariable is the number
        %   of variables in the model.
        %
        %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the network model at
        %   X sampled points, where X can be of size greater than nVariable.
        %   VAROBJ will specify which columns of X to be used in evaluating
        %   the model to produce a column vector Y of the model's output.
        
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         xT = obj.xTransform;
         yT = obj.yTransform;
         weights = obj.Weights;
         biases = obj.Biases;
         actfcn = obj.ActivationFcn;
         nS = size(X,1);
         X = 2*(X-repmat(xT.Translation,nS,1))./repmat(xT.Range,nS,1)-1;
         if size(weights{1},2) ~= size(X,2)
            error('Wrong input dimension of variables')
         else
            actf = actfcn{1};
            y = actf(X*weights{1}'+repmat(biases{1}',nS,1));
            for i = 2:obj.LayerNumber
               actf = actfcn{i};
               y = actf(y*weights{i}'+repmat(biases{i}',nS,1));
            end
         end
         y = 0.5*(y+1)*yT.Range+yT.Translation+obj.Constant;
      end
      
      function y = get.LayerNumber(obj)
         y = length(obj.Weights);
      end

   end
   
      
   
   methods (Hidden = true)
      function f = calculateGradient(obj)
         f = @gradFcn;
         nLayer = obj.LayerNumber;
         weights = obj.Weights;
         biases = obj.Biases;
         actfcn = obj.ActivationFcn;
         gfcn = obj.ActivationGradientFcn;
         yRange = obj.yTransform.Range;
         xRange = obj.xTransform.Range;
         xTrans = obj.xTransform.Translation;       
         function g = gradFcn(x)
            z = cell(nLayer,1);
            xx = 2*(x-xTrans')./xRange'-1;
            for i = 1:nLayer
               ww = weights{i};
               z{i} = ww*xx+biases{i};
               xx = actfcn{i}(z{i});
            end
            gg = 0.5*yRange;
            for i = nLayer-1:-1:1
               ww = weights{i+1};
               n1 = size(ww,2);
               g = zeros(n1,1);
               for j = 1:length(gg)
                  g = g+gg(j)*ww(j,:)'.*gfcn{i}(z{i});
               end
               gg = g;
            end
            ww = weights{1};
            n1 = size(ww,2);
            g = zeros(n1,1);
            for j = 1:length(gg)
               g = g + gg(j)*ww(j,:)';
            end
            g = 2*g./xRange';
         end
      end
   end
   
end

