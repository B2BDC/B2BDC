classdef VariableList < B2BDC.Util.IContainer
% A container for B2BDC.B2Bvariables.ModelVariable objects
   
% Created: June 17, 2015, myf + Wenyu Li
% Modified: Nov 1, 2015     Wenyu Li (Sample function added)
   
   
   properties(Dependent, Hidden = true)
      TransMatrix = []; % Transform matrix to transform x from unnormalized space specified by Variables to normalized space ([-1,1])
   end
   
   properties (SetAccess = private)
      ExtraLinConstraint = struct('A',[],'LB',[],'UB',[]);
<<<<<<< Updated upstream
=======
      ExtraQuaConstraint = struct('Q',[],'UB',[],'xStart',[]);
>>>>>>> Stashed changes
      ScalingVector = [];
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
         
         A1 = obj.ExtraLinConstraint.A;
         LB1 = obj.ExtraLinConstraint.LB;
         UB1 = obj.ExtraLinConstraint.UB;
         if isempty(obj.Values)
            name1 = {};
         else
            name1 = {obj.Values.Name};
         end
         A2 = varList.ExtraLinConstraint.A;
         LB2 = varList.ExtraLinConstraint.LB;
         UB2 = varList.ExtraLinConstraint.UB;
         name2 = {varList.Values.Name};
         obj = obj.clearExtraConstraint;
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
         allName = {obj.Values.Name};
         nV = obj.Length;
         if ~isempty(A1)
            n1 = size(A1,1);
            At = zeros(n1,nV);
            [~,~,id] = intersect(name1,allName,'stable');
            At(:,id) = A1;
            obj = obj.addLinearConstraint(At,LB1,UB1);
%             obj = obj.addLinearConstraint([At;-At],[UB1;-LB1]);
         end
         if ~isempty(A2)
            n2 = size(A2,1);
            At = zeros(n2,nV);
            [~,~,id] = intersect(name2,allName,'stable');
            At(:,id) = A2;
            obj = obj.addLinearConstraint(At,LB2,UB2);
%             obj = obj.addLinearConstraint([At;-At],[UB2;-LB2]);
<<<<<<< Updated upstream
         end 
=======
         end
         if ~obj.checkFeasibility
            disp('The resulted variableList is not feasible')
         end
>>>>>>> Stashed changes
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
         A = obj.ExtraLinConstraint.A;
         tA = A(:,id);
         if any(tA)
            error('Extra linear constraints contain the variables you want to delete')
         end
         if ~isempty(id)
            obj = obj.remove(id);
            obj.ScalingVector = [];
         else
            disp('The variable you want to delete is not in the current list')
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
            if ~isempty(obj.ExtraLinConstraint.A)
               A = obj.ExtraLinConstraint.A;
               LB = obj.ExtraLinConstraint.LB;
               UB = obj.ExtraLinConstraint.UB;
               Var = Var.addLinearConstraint(A,LB,UB);
%                Var = Var.addLinearConstraint([A;-A],[UB;-LB]);
            end
            obj.ScalingVector = [];
         else
            error('Wrong number of input arguments')
         end
      end
      
      function newVar = addLinearConstraint(obj,A,b,c)
         % NEWVAR = ADDEXTRALINCONSTRAINT(OBJ,A,B,C) returns the new
         % VariableList object with added linear constraints on the
         % variables satisfying A*x <= b if nargin < 4. Otherwise returns
         % the new VariableList object with added linear constraints that
         % b <= A*x <= c.
         
         if nargin < 4
            if size(b,1) == 1
               b = b';
            end
            newVar = obj;
            [nC, nVar] = size(A);
            if nVar ~= obj.Length
               error('Invalid input matrix dimension')
            end
            if nC ~= length(b)
               error('Invalid input vector dimension')
            end
            [idA,id] = intersect(A,-A,'rows','stable');
            nint = 0.5*size(idA,1);
            nuni = nC - 2*nint;
            A0 = zeros(nint+nuni,nVar);
            LB = zeros(nint+nuni,1);
            UB = zeros(nint+nuni,1);
            if ~isempty(idA)
               count = 1;
               for i = 1:2*nint
                  idx = find(idA(i,:),1);
                  if idA(i,idx) > 0
                     A0(count,:) = idA(i,:);
                     [~,id1] = intersect(A,idA(i,:),'rows');
                     [~,id2] = intersect(-A,idA(i,:),'rows');
                     UB(count) = b(id1);
                     LB(count) = -b(id2);
                     count = count+1;
                  end
               end
            end
            tmpA = A;
            tmpb = b;
            A(id,:) = [];
            b(id,:) = [];
            A0(nint+1:end,:) = A;
            UB(nint+1:end) = b;
            for i = 1:nuni
               tmpC = A(i,:);
               LB(nint+i) = finddirectionmin(obj,tmpA,tmpb,tmpC);
            end
            for i = 1:nuni
               idx = find(A0(nint+i,:),1);
               if A0(nint+i,idx) < 0
                  A0(nint+i,:) = -A(i,:);
                  tmpLB = -UB(nint+i);
                  UB(nint+i) = -LB(nint+i);
                  LB(nint+i) = tmpLB;
               end
            end
            if ~isempty(obj.ExtraLinConstraint.A)
               Ai = obj.ExtraLinConstraint.A;
               Li = obj.ExtraLinConstraint.LB;
               Ui = obj.ExtraLinConstraint.UB;
               [~,id1,id2] = intersect(Ai,A0,'rows','stable');
               if ~isempty(id1)
                  %                LB(LB(id2) < Li(id1)) = Li(id1);
                  LB(id2) = max([LB(id2), Li(id1)],[],2);
                  %                UB(UB(id2) > Ui(id1)) = Ui(id1);
                  UB(id2) = min([UB(id2), Ui(id1)],[],2);
               end
               Ai(id1,:) = [];
               Li(id1) = [];
               Ui(id1) = [];
               A0 = [A0;Ai];
               LB = [LB;Li];
               UB = [UB;Ui];
            end
            newVar.ExtraLinConstraint.A = A0;
            newVar.ExtraLinConstraint.LB = LB;
            newVar.ExtraLinConstraint.UB = UB;
         else
            if size(A,1) ~= length(b)
               error('Wrong lower bound dimension')
            end
            if size(A,1) ~= length(c)
               error('Wrong upper bound dimension')
            end
            newVar = obj;
<<<<<<< Updated upstream
            newVar.ExtraLinConstraint.A = A;
            newVar.ExtraLinConstraint.LB = b;
            newVar.ExtraLinConstraint.UB = c;
         end
         newVar.ScalingVector = [];
=======
            A0 = obj.ExtraLinConstraint.A;
            b0 = obj.ExtraLinConstraint.LB;
            c0 = obj.ExtraLinConstraint.UB;
            newVar.ExtraLinConstraint.A = [A;A0];
            newVar.ExtraLinConstraint.LB = [b;b0];
            newVar.ExtraLinConstraint.UB = [c;c0];
         end
         newVar.ScalingVector = [];
         if ~newVar.checkFeasibility
            disp('The resulted variableList is not feasible')
            newVar = obj;
         end
>>>>>>> Stashed changes
      end
      
      function newVar = makeSubset(obj,varIdx)
         % NEWVAR = MAKESUBSET(OBJ,INDEX) returns a VariableList object
         % including only a subset of the original VariableList object
         % specifed by the index input
         
         vName = {obj.Values.Name}';
         if ischar(varIdx)
            [~,~,id] = intersect(varIdx,{obj.Values.Name},'stable');
         elseif isvector(varIdx)
            id = varIdx;
         else
            error('Wrong input index type or length')
         end
<<<<<<< Updated upstream
         if ~isempty(obj.ExtraLinConstraint.A)
            A = obj.ExtraLinConstraint.A;
            tA = A(:,id);
            if any(tA)
               error('Extra linear constraints contain variables outside the subset')
            end
         end
=======
>>>>>>> Stashed changes
         H = obj.calBound;
         ob = [obj.Values.NominalValue]';
         newName = vName(id);
         newH = H(id,:);
         newOB = ob(id,:);
         newVar = generateVar(newName, newH, newOB);
         if ~isempty(obj.ExtraLinConstraint.A)
            A = obj.ExtraLinConstraint.A;
            LB = obj.ExtraLinConstraint.LB;
            UB = obj.ExtraLinConstraint.UB;
<<<<<<< Updated upstream
            A = A(:,id);
            idA = true(size(A,1),1);
            for i = 1:size(A,1)
               tmpid = find(A(i,:));  
               if length(tmpid) > 1
                  idA(i) = false;
               end
            end
            A(idA,:) = [];
            if ~isempty(A)
               LB(idA) = [];
               UB(idA) = [];
               newVar = newVar.addLinearConstraint(A,LB,UB);
=======
            A_new = A(:,id);
            idA = true(size(A,1),1);
            for i = 1:size(A,1)
               a1 = sum(A(i,:)~=0);
               a2 = sum(A_new(i,:)~=0);
               if a1 == a2
                  idA(i) = false;
               end
            end
            A_new(idA,:) = [];
            if ~isempty(A)
               LB(idA) = [];
               UB(idA) = [];
               newVar = newVar.addLinearConstraint(A_new,LB,UB);
>>>>>>> Stashed changes
%                newVar = newVar.addLinearConstraint([A;-A],[UB;-LB]);
            end
         end
      end
      
      function newVar = clearExtraConstraint(obj)
<<<<<<< Updated upstream
         newVar = obj;
         newVar.ExtraLinConstraint = struct('A',[],'LB',[],'UB',[]);
         newVar.ScalingVector = [];
      end
      
=======
         % clear extra constraints
         newVar = obj;
         newVar.ExtraLinConstraint = struct('A',[],'LB',[],'UB',[]);
         newVar.ExtraQuaConstraint = struct('Q',[],'UB',[],'xStart',[]);
         newVar.ScalingVector = [];
      end
      
      function y = checkFeasibility(obj) 
         y = true;
         if isempty(obj.ExtraQuaConstraint.Q)
            if ~isempty(obj.ExtraLinConstraint.A)
               warning('off','all');
               y = false;
               for i = 1:10
                  opt1 = optimoptions('linprog');
                  opt1.Display = 'none';
                  A = obj.ExtraLinConstraint.A;
                  ub = obj.ExtraLinConstraint.UB;
                  lb = obj.ExtraLinConstraint.LB;
                  H = obj.calBound;
                  x0 = linprog(zeros(size(A,2),1),[A;-A],[ub;-lb],[],[],H(:,1),...
                     H(:,2),opt1);
                  if ~isempty(x0)
                     y = true;
                     break
                  end
               end
            end
         else
            x0 = obj.ExtraQuaConstraint.xStart;
            if ~obj.isFeasiblePoint(x0')
               y = false;
            end
         end
      end
      
      function newVar = addQuadraticConstraint(obj,Q,UB,x0)
         newVar = obj;
         if size(x0,1) == 1
            x0 = x0';
         end
         if ~obj.isFeasiblePoint(x0')
            msgbox('The provided point is not feasible for the original varList')
            return
         end
         Q = 0.5*(Q+Q');
         [~,d] = eig(Q);
         if any(diag(d) <= -1e-14)
            msgbox('The provided quadratic constraint is not convex')
            return
         end
         if [1;x0]'*Q*[1;x0] <= UB
            newVar.ExtraQuaConstraint.Q{end+1} = Q;
            newVar.ExtraQuaConstraint.UB(end+1) = UB;
            newVar.ExtraQuaConstraint.xStart = x0;
         end
      end
      
   end
   
   methods (Static, Hidden = true)
       xval = q2sample(Mgd,idx,H,Xvals,V);
       xVal = quadSample(Q,H,x0,V)
>>>>>>> Stashed changes
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
         if ~isempty(obj.ExtraLinConstraint.A)
            obj.ExtraLinConstraint.A = [A, zeros(size(A,1),1)];
         end
      end
   end
   
   
end