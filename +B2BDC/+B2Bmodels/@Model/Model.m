classdef Model < handle
    % An abstract super class for B2BDC models
    
    % Created: June 17, 2015   Wenyu Li
    % Modified: Nov 1, 2015     Wenyu Li (error stats added)
    
    properties
        Variables %A VariableList object contains all variables associated with this model
    end
    
    properties
        ErrorStats = struct('absMax',[],'absAvg',[],'relMax',[],'relAvg',[]); % Statistical information about the fitting error (if any) 
    end
    
    properties (Dependent)
        VarRanges % A nVar-by-2 matrix represents the range of all variables associated with the model
        VarNames  % A cell array including the nameof all variables associated with the model 
    end
    
    methods
        function y = get.VarRanges(obj)
            %   Y = VARRANGES(OBJ) returns a matrix Y of size nVariable-by-2
            %   containing the range of all variables in the model. 
            
            modelvar = obj.Variables.Values;
            n_variable = length(modelvar);
            y = zeros(n_variable,2);
            for i = 1:n_variable
                variable = modelvar(i);
                y(i,1) = variable.LowerBound;
                y(i,2) = variable.UpperBound;
            end
        end
        
        function y = get.VarNames(obj)
            %   Y = VARNAMES(OBJ) returns a cell array of size
            %   nVariable-by-1 containing the name of all variables in the
            %   model.
            
            modelvar = obj.Variables.Values;
            n_variable = length(modelvar);
            y = cell(n_variable,1);
            for i = 1:n_variable
                variable = modelvar(i);
                y{i} = variable.Name;
            end
        end
        
    end
    
    methods (Abstract)
        %   Y = EVAL(OBJ, X) evaulates a model at X sampled points to
        %   produce a column vector Y of the model's output.
        
        y = eval(obj,X,varObj)
    end
    
    methods (Hidden = true)
       function obj = addVar(obj,variable)
            if isempty(obj.Variables)
                obj.Variables = B2BDC.B2Bvariables.VariableList();
            end
            obj.Variables = obj.Variables.add(variable);
        end
    end
    
    methods (Hidden = true)
      % Matlab built-in functions
       y = addlistener()
       y = findobj() 
       y = findprop() 
       y = notify() 
       y = le()
       y = lt() 
       y = ne()
       y = ge() 
       y = gt() 
       y = eq() 
   end
    
end