classdef ModelVariable
   % A ModelVariable object contains the information of one variable
   
   % Created: June 17, 2015
   
   properties
      NominalValue  % Nominal value of the variable
      LowerBound    % Lower bound of the variable
      UpperBound    % Upper bound of the variable
      Name          % Descriptive name of the variable
   end
   
   methods
      function obj = ModelVariable(name,LB,UB,val)
         %   OBJ = MODELVARIABLE(NAME, LB, UB) returns a 
         %   B2BDC.B2Bvariables.ModelVariable object. The name is a
         %   character string describing the variable, LB and UB are the
         %   lower and upper bounds to the variable. The NominalValue will 
         %   taken as the mean of the lower and upper bounds. 
         % 
         %   OBJ = MODELVARIABLE(NAME, LB, UB, VAL) the NominalValue of the 
         %   ModelVariable can be specified by VAL.
         
         if nargin > 0
            if ischar(name)
               obj.Name = name;
            else
               error('Variable name should be a string');
            end
            if ~isscalar(UB) || ~isscalar(LB)
               error('Uncertainty bound should be scalar')
            elseif LB >= UB
               error(['The domain of variable ' name ' is empty']) 
            else
                obj.LowerBound = LB;
                obj.UpperBound = UB;
            end
            if nargin > 3
               obj.NominalValue = val;
            else
               obj.NominalValue = 0.5*(LB+UB);
            end
         end
      end
   end
   
end

