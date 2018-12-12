classdef OptimOption < handle
   % A B2BDC.Option.OptimOption object that specifies the optimization
   % criteria for finding optimal parameters. The syntax of the input and 
   % the meaning of each property are in the following:
   % opt = B2BDC.Option.OptimOption(PropName, PropValue ... )
   % ---------------------------------------------------------------------
   % OptimizationMethod:
   %   LSF(default) - Least square optimization in feasible set
   %   LSH - Least square optimization in prior H
   %   1NF - 1 norm optimization in F
   % ---------------------------------------------------------------------
   % PenaltyWeight:
   %   relative(default) - squared relative difference
   %   absolute - squared absolute difference
   %   user-defined - user defined weights
   % ---------------------------------------------------------------------
   % RandomStart:
   %   A positive integer specifying the number of starting points, the
   %   default value is 1.
   % ---------------------------------------------------------------------
   
   % Created: Oct 25, 2017    Wenyu Li
   
   
   properties
      OptimizationMethod = [];
      PenaltyWeight = [];
      RandomStart = [];
      PredictionTol = [];
      VariableTol = [];
   end
   
   methods
      function  obj = OptimOption(inputCell)
         % To generate a B2BDC.Option object
         p = {'OptimizationMethod','PenaltyWeight','RandomStart'};
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
                        obj.OptimizationMethod = inputCell{2*i};
                     case 2
                        obj.PenaltyWeight = inputCell{2*i};
                     case 3
                        obj.RandomStart = inputCell{2*i};
                     case 4
                        obj.PredictionTol = inputCell{2*i};
                     case 5
                        obj.VariableTol = inputCell{2*i};
                  end
               else
                  error('Invalid input property names')
               end
            end
         end
         if isempty(obj.OptimizationMethod)
            obj.OptimizationMethod = 'LSF';
         end
         if isempty(obj.PenaltyWeight)
            obj.PenaltyWeight = 'relative';
         end
         if isempty(obj.RandomStart)
            obj.RandomStart = 10;
         end
         if isempty(obj.RandomStart)
            obj.RandomStart = 1;
         end
         if isempty(obj.PredictionTol)
            obj.PredictionTol = 0;
         end
         if isempty(obj.VariableTol)
            obj.VariableTol = 0;
         end
      end
      
      function set.OptimizationMethod(obj,str1)
         c = {'LSF','LSH','1NF'};
         if any(strcmp(c,str1))
            obj.OptimizationMethod = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.PenaltyWeight(obj,str1)
         c = {'relative','absolute','user-defined'};
         if any(strcmp(c,str1))
            obj.PenaltyWeight = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.RandomStart(obj,c1)
         if round(c1) < 1
            obj.RandomStart = 1;
         else
            obj.RandomStart = round(c1);
         end
      end
      
      function set.VariableTol(obj,c)
         if c>=0 && c<1
            obj.VariableTol = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.PredictionTol(obj,c)
         if c>=0 && c<1
            obj.PredictionTol = c;
         else
            error('Invalid input property value')
         end
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