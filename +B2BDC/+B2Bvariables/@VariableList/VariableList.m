classdef VariableList < B2BDC.Util.IContainer
% A container for B2BDC.B2Bvariables.ModelVariable objects
   
% Created: June 17, 2015, myf + Wenyu Li
% Modified: Nov 1, 2015     Wenyu Li (Sample function added)
   
   
   properties(Dependent, Hidden = true)
      TransMatrix = []; % Transform matrix to transform x from unnormalized space specified by Variables to normalized space ([-1,1])
   end
   
   methods
      
      function obj = VariableList()
         %   OBJ = VARIABLELIST() returns a B2BDC.B2Bvariables.VariableList
         %   object, which acts as a container for the ModelVariable 
         %   objects.
         
         valClass = 'B2BDC.B2Bvariables.ModelVariable';
         obj = obj@B2BDC.Util.IContainer(valClass);
      end
      
      function obj = addList(obj,varList)
         %   OBJ = ADDLIST(OBJ, VARLIST) will return a VariableList OBJ 
         %   which has been added to another VariableList VARLIST. If any 
         %   Variables are in common between the two lists (based on the 
         %   Variable Name) the intersection of their LowerBound and 
         %   UpperBound will be taken.
         
         variables = varList.Values;
         for i = 1:varList.Length
            variable = variables(i);
            [y,ind] = obj.find('Name',variable.Name);
            if isempty(y)
               obj = obj.add(variable);
            else
               if variable.LowerBound > obj.Values(ind).LowerBound
                  obj.Container{ind}.LowerBound = variable.LowerBound;
               end
               if variable.UpperBound < obj.Values(ind).UpperBound
                  obj.Container{ind}.UpperBound = variable.UpperBound;
               end
            end
         end
      end
      
      function xSample = makeLHSsample(obj,nSample)
         %   XSAMPLE = MAKELHSSAMPLE(OBJ, NSAMPLE) returns latin hypercube 
         %   samples of the VariableList OBJ, within the Variable 
         %   LowerBound and UpperBound. NSAMPLE specifies the number of 
         %   sample points to return from the VariableList domain.
         
         nVar = obj.Length;
         H = [[obj.Values.LowerBound]', [obj.Values.UpperBound]'];
         dH = diff(H');
         xdesign = lhsdesign(nSample,nVar);
         xSample = repmat(H(:,1)',nSample,1) + repmat(dH,nSample,1).*xdesign;
      end
      
      function y = get.TransMatrix(obj)
         nVar = obj.Length;
         mx = mean(obj.calBound')';
         dx = 0.5*diff(obj.calBound')';
         y = [1, zeros(1,nVar);
            -mx./dx, diag(1./dx)];
      end

      function obj = deleteVariable(obj,varIdx)
          %   OBJ = DELETEVARIABLE(OBJ, VARIDX) removes specified 
          %   variable from the VariableList and returns a VariableList
          %   OBJ. VARIDX will specify which variable to delete by either 
          %   a character string of the Variable name or the index of the
          %   Variable in the VariableList.
          
         if ischar(varIdx)
            [~,id,~] = intersect(varIdx,{obj.Values.Name});
         elseif isscalar(varIdx)
            id = varIdx;
         else
            error('Wrong input index type or length')
         end
         if ~isempty(id)
            obj = obj.remove(id);
         else
            disp('The variable you want to delete is not in the current list')
         end
      end
      
      function obj = scale(obj,factor)
         %   OBJ = SCALE(OBJ, FACTOR) returns a VariableList OBJ where all 
         %   variables have been scaled by FACTOR. FACTOR can be a
         %   vector of length two where the lower bound is scaled by the
         %   first element and the upper bound is scaled by the scond.
         %   FACTOR also can be a nVar-by-2 matrix which scales each 
         %   Variable LowerBound and UpperBound by the corresponding 
         %   elements of the FACTOR. 
         
         if isvector(factor)
            if length(factor) ~= 2
               error(['Incorrect dimensions for a scaling factor. ' ...
                   'Scaling factor should be a vector of length two or ' ...
                   'an nVar-by-2 matrix.'])
            else
               for i = 1:obj.Length
                  obj.Container{i}.LowerBound = factor(1)*obj.Container{i}.LowerBound;
                  obj.Container{i}.UpperBound = factor(2)*obj.Container{i}.UpperBound;
               end
            end
         elseif ismatrix(factor)
            if size(factor,1) ~= obj.Length || size(factor,2) ~= 2
               error(['Incorrect dimensions for a scaling factor. ' ...
                   'Scaling factor should be a vector of length two or ' ...
                   'an nVar-by-2 matrix.'])
            else
               for i = 1:obj.Length
                  obj.Container{i}.LowerBound = factor(i,1)*obj.Container{i}.LowerBound;
                  obj.Container{i}.UpperBound = factor(i,2)*obj.Container{i}.UpperBound;
               end
            end
         else
            error(['Incorrect input for a scaling factor. ' ...
                 'Scaling factor should be a vector of length two or ' ...
                 'an nVar-by-2 matrix.'])
         end
         bds = obj.calBound;
         if any(bds(:,1) >= bds(:,2))
            error('The domain after scaling cannot be empty')
         end
      end
      
      function varBD = calBound(obj)
         %   VARBD = CALBOUND(OBJ) returns a nVar-by-2 matrix containing 
         %   the LowerBound and UpperBound of all Variables (nVar) in the
         %   VariablesList OBJ.
         
         varVal = obj.Values;
         varBD = [[varVal.LowerBound]', [varVal.UpperBound]'];
      end
  
      function nVar = length(obj)
         %   NVAR = LENGTH(OBJ) returns the number of Variables in the 
         %   VariableList OBJ. 
         
          nVar = obj.Length;
      end
      
      function nVar = numel(obj)
         %   NVAR = NUMEL(OBJ) returns the number of Variables in the 
         %   VariableList OBJ. 
         
          nVar = obj.Length;
      end
      
      function flag = isSubset(obj, var2)
         % Returns true if the domain specified by obj is a subset of the
         % domain specified by var2, otherwise false.
         bd1 = obj.calBound;
         bd2 = var2.calBound;
         flag = true;
         if any(bd1(:,1) < bd2(:,1))
            flag = false;
         elseif any(bd1(:,2) > bd2(:,2))
            flag = false;
         end   
      end
      
      function domOut = splitToinclude(obj,var2)
         % Calculate the complementary of domain specified by var2. The
         % whole domain is specified by obj.
         if var2.isSubset(obj)
            varNames = {obj.Values.Name};
            tgtH = var2.calBound;
            domH = obj.calBound;
            domOut = {};
            for i = 1:length(varNames)
               if tgtH(i,1) > domH(i,1)
                  tmpH = domH;
                  tmpH(i,2) = tgtH(i,1);
                  domH(i,1) = tgtH(i,1);
                  domOut{end+1} = generateVar(varNames, tmpH);
               end
               if tgtH(i,2) < domH(i,2)
                  tmpH = domH;
                  tmpH(i,1) = tgtH(i,2);
                  domH(i,2) = tgtH(i,2);
                  domOut{end+1} = generateVar(varNames, tmpH);
               end
            end
            domOut = [domOut{:}]';
         else
            error('The target domain is not fully contained in the whole domain')
         end
      end
      
      function varInt = findIntersect(obj, var2)
         % Returns a VariableList that specifies the domain as the
         % intersection of the two domains specified by obj and var2
         % respectively. It returns empty if the two domains don't
         % intersect.
         varNames = {obj.Values.Name};
         h1 = obj.calBound;
         h2 = var2.calBound;
         id1 = h1(:,1) >= h2(:,1);
         h1(~id1,1) = 0;
         h2(id1,1) = 0;
         id2 = h1(:,2) <= h2(:,2);
         h1(~id2,2) = 0;
         h2(id2,2) = 0;
         hInt = h1 + h2;
         if any(hInt(:,1) >= hInt(:,2))
            varInt = [];
         else
            varInt = generateVar(varNames, hInt);
         end
      end
      
      function Var = changeBound(obj, newBD, idx)
         % Modifies the bounds of specified VariableList
         % The input arguments are:
         %  newBD - A nChange-by-2 (or 3) matrix defines the new
         %          bounds (new observed values)
         %  idx - A cell array of variable names or an index array of
         %        length nChange that specifies the corresponding variables
         %        with repsect to the new bounds. If this input is
         %        not given and the newBD matches the length of the
         %        VariableList, it modifies all variables in the
         %        VariableList
         allVarName = {obj.Values.Name};
         if nargin == 3
            if length(idx) ~= size(newBD,1)
               error('Inconsistent new bounds size with number of variables')
            elseif iscell(idx)
               for i = 1:length(idx)
                  [~,id] = intersect(allVarName, idx{i});
                  if isempty(id)
                     error(['The variable ' idx{i} ' is not in the VariableList'])
                  end
                  varName = obj.Values(id).Name;
                  if size(newBD,2) == 2
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2));
                  else
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2), newBD(i,3));
                  end
                  Var = obj.replace('Name',varName,newVar);
               end
            elseif isnumeric(idx)
               for i = 1:length(idx)
                  varName = obj.Values(idx(i)).Name;
                  if size(newBD,2) == 2
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2));
                  else
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2), newBD(i,3));
                  end
                  Var = obj.replace('Name',varName,newVar);
               end
            else
               error('Wrong input index')
            end
         elseif nargin == 2 && size(newBD,1) == obj.Length
            varName = {obj.Values.Name};
            if size(newBD,2) == 2
               Var = generateVar(varName,newBD);
            else
               Var = generateVar(varName,newBD(:,1:2),newBD(:,3));
            end
         else
            error('Wrong number of input arguments')
         end
      end
 
   end
   
   methods (Hidden = true)
      function obj = addVariable(obj,varName,varRange,varValue)
         if nargin > 3
            variable = B2BDC.B2Bvariables.ModelVariable(varName,varRange(1),varRange(2),varValue);
         else
            variable = B2BDC.B2Bvariables.ModelVariable(varName,varRange(1),varRange(2));
         end
         [y,ind] = obj.find('Name',variable.Name);
         if isempty(y)
            obj = obj.add(variable);
         else
            if variable.LowerBound > obj.Values(ind).LowerBound
               obj.Container{ind}.LowerBound = variable.LowerBound;
            end
            if variable.UpperBound < obj.Values(ind).UpperBound
               obj.Container{ind}.UpperBound = variable.UpperBound;
            end
         end
      end
   end
   
end