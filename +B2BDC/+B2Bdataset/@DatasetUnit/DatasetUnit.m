classdef DatasetUnit
    % A DatasetUnit belonging to the B2BDC Dataset
   
   % Created: June 17, 2015     myf,  Wenyu Li
   % Modified: December 11, 2015    Jim Oreluk - Added error message if
   % ub < lb
   
   properties
      Name   % Descriptive name for the dataset unit
        LowerBound   % Lower bound of the observed experimental value
        ObservedValue  % Experimental Value
        UpperBound   % Upper bound of the observed experimental value
        SurrogateModel  % Surrogate model representing the experiment
   end
   
   properties (Dependent)
      VariableList  % VariableList object containing all model variables associated with this dataset unit
   end
   
   methods
      function obj = DatasetUnit(dsUnitname,modelObj,LB,UB,val)
         % Constructor of the B2BDC.B2Bdataset.DatasetUnit object.
            %
            % The input arguments are:
            %   name - a string defining a descriptive name for the dataset unit
            %   modelObj - A B2BDC.B2Bmodels.Model subclass
            %   LowerBound - Lower bound on the ObservedValue
            %   UpperBound - Upper bound on the ObservedValue
            %   ObservedValue - Experimental observation (optional). If not defined, ObservedValue is calculated as the mean of LowerBound and UpperBound.
         
         if nargin > 0
            obj.Name = dsUnitname;
            if ~isa(modelObj,'B2BDC.B2Bmodels.Model')
               error('Model must be a Model class object.')
            else
               obj.SurrogateModel = modelObj;
            end
            if ~isscalar(UB) || ~isscalar(LB)
               error('Upper and lower bounds should be scalar values.')
            else
               obj.LowerBound = LB;
               obj.UpperBound = UB;
            end
            if nargin > 4
               obj.ObservedValue = val;
            else
               obj.ObservedValue = 0.5*(LB+UB);
            end
         end
      end
      
      function y = get.VariableList(obj)
         %   Y = VARIABLELIST(OBJ) returns the B2BDC.B2Bvariable.VariableList
         %   from the OBJ.SurrogateModel.
         
         model = obj.SurrogateModel;
         y = model.Variables;
      end
      
      function y = changeBounds(obj,newbd)
          %   Y = CHANGEBOUNDS(OBJ, NEWBD) changes the lower and upper bounds
          %   of OBJ by the vector NEWBD. When NEWBD is of length 2, OBJ
          %   will be modified as [OBJ.LowerBound, OBJ.UpperBound].
          %   OBJ.ObservedValue will be changed to the mean value the
          %   LowerBound and UpperBound.
          %
          %   When NEWDB is of length 3, OBJ will be modified as
          %   [OBJ.LowerBound, OBJ.UpperBound, OBJ.ObservedValue].
          
         lb = newbd(1);
         ub = newbd(2);
         if ub < lb
             error('Upper bound cannot be less than the lower bound.')
         end
         dsName= obj.Name;
         dsMod = obj.SurrogateModel;
         if length(newbd) < 3
            ob = 0.5*(lb+ub);
         else
            ob = newbd(3);
         end
         y = B2BDC.B2Bdataset.DatasetUnit(dsName,dsMod,lb,ub,ob);
      end
   end
   
   methods (Hidden = true)
      flag = quadratictest(dsUnits)
      y = lintest(dsUnits)
   end
   
end

