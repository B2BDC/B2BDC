classdef SOSOption < handle
   % A B2BDC.Option.POPOption object that specifies the SparsePOP criteria.
   % The syntax of the input and the meaning of each property are in the
   % following:
   % opt = B2BDC.Option.SampleOption(PropName, PropValue ... )
   % More detail information can be found in SparsePOP toolbox.
   
   
   % Created: Oct 11, 2016    Wenyu Li
   %  Modified: Septempber 2, 2015   Wenyu Li (Data normalization option added)
   
   
   properties
      SOSOrder = [];
      sparseSW = [];
      perturbation = [];
      SDPsolver = [];
      printFileName = [];
      printLevel = [];
      POPsolver = [];
      mex = [];
   end
   
   properties (Dependent)
      Value = [];
   end
   
   methods
      function  obj = POPOption(inputCell)
         % To generate a B2BDC.Option object
         p = {'SOSOrder','sparseSW','perturbation',...
            'SDPsolver','printFileName','printLevel',...
            'POPsolver','mex'};
         if nargin > 0 
            nin = length(inputCell);
         else
            nin = 0;
         end
         if mod(nin,2) ~= 0
            error('Wrong number of input argument')
         end
         nset = floor(0.5*nin);
         if nset > 0
            for i = 1:nset
               if any(strcmp(p,inputCell{2*i-1}))
                  id = find(strcmp(p,inputCell{2*i-1}) == 1);
                  switch id
                     case 1
                        obj.SOSOrder = inputCell{2*i};
                     case 2
                        obj.sparseSW = inputCell{2*i};
                     case 3
                        obj.perturbation = inputCell{2*i};
                     case 4
                        obj.SDPsolver = inputCell{2*i};
                     case 5
                        obj.printFileName = inputCell{2*i};
                     case 6
                        obj.printLevel = inputCell{2*i};
                     case 7
                        obj.POPsolver = inputCell{2*i};
                     case 8
                        obj.mex = inputCell{2*i};
                  end
               else
                  error('Invalid input property names')
               end
            end
         end
         if isempty(obj.relaxOrder)
            obj.SOSOrder = 0;
         end
         if isempty(obj.sparseSW)
            obj.sparseSW = 1;
         end
         if isempty(obj.perturbation)
            obj.perturbation = 0;
         end
         if isempty(obj.SDPsolver)
            obj.SDPsolver = 'sedumi';
         end
         if isempty(obj.printFileName)
            obj.printFileName = 0;
         end
         if isempty(obj.printLevel)
            obj.printLevel = [0,0];
         end
         if isempty(obj.POPsolver)
            obj.POPsolver = 'interior-point';
         end
         if isempty(obj.mex)
            obj.mex = 1;
         end
      end
      
      function set.SOSOrder(obj,num1)
         if num1 >= 0
            obj.SOSOrder = ceil(num1);
         else
            error('Invalid input property value')
         end
      end
      
      function set.sparseSW(obj,num1)
         if num1 == 0 || num1 == 1
            obj.sparseSW = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.perturbation(obj,num1)
         if num1 >= 0
            obj.perturbation = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.SDPsolver(obj,str1)
         c = {'sedumi','sdpa'};
         if any(strcmp(c,str1))
            obj.SDPsolver = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.printFileName(obj,num1)
         if num1 == 0 || num1 == 1
            obj.printFileName = num1;
         elseif ischar(num1)
            obj.printFileName = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.printLevel(obj,num1)
         if length(num1) == 1
            if num1 == 0 || num1 == 1 || num1 == 2
               obj.printLevel(1) = num1;
            else
               error('Invalid input property value')
            end
         elseif length(num1) == 2
            c1 = num1(1);
            c2 = num1(2);
            if c1 == 0 || c1 == 1 || c1 == 2
               obj.printLevel(1) = c1;
            else
               error('Invalid input property value')
            end
            if c2 == 0 || c2 == 1 || c2 == 2
               obj.printLevel(2) = c2;
            else
               error('Invalid input property value')
            end
         else
            error('Invalid input property value')
         end
      end
      
      function set.POPsolver(obj,str1)
         c = {'active-set','trust-region-reflective',...
            'interior-point'};
         if any(strcmp(c,str1))
            obj.POPsolver = str1;
         else
            obj.POPsolver = [];
         end
      end
      
      function set.mex(obj,num1)
         if num1 == 0 || num1 == 1
            obj.mex = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function y = get.Value(obj)
         y.SOSOrder = obj.SOSOrder;
         y.sparseSW = obj.sparseSW;
         y.perturbation = obj.perturbation;
         y.SDPsolver = obj.SDPsolver;
         y.printFileName = obj.printFileName;
         y.printLevel = obj.printLevel;
         y.POPsolver = obj.POPsolver;
         y.mex = obj.mex;
      end
   end
   
   methods (Hidden = true)
      % Matlab additional Functions
       y = addlistener() %matlab add listener events
       y = delete() % matlab delete handle
       y = findobj() % matlab find objects with a specified property value
       y = findprop() % matlab find matlab property of handle obj
       y = notify() 
       y = le() % matlab less than or equal to
       y = lt() % matlab less than
       y = ne() % matlab not equal to
       y = ge() % matlab greater than or equal to
       y = gt() % matlab greater than
       y = eq() % matlab equal to
%        y = isvalid() % matlab method for timer obj -- breaks when
       %hidden...
   end
   
end