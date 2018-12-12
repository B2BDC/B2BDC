classdef Dataset < handle
   % B2BDC Dataset Object
   % A dataset is a collection of models, observations, variables and their
   % respective bounds. 

   %  Created: June 19, 2015   Wenyu Li
   %  Modified: June 22, 2015    Wenyu Li
   %  Modified: July 2, 2015    Wenyu Li  (option added for consistency
   %                                       measure)
   %  Modified: July 5, 2015    Wenyu Li  (Sensitivity and Range property
   %                                       for QOI added)
   %  Modified: Oct 22, 2015    Wenyu Li  (feasible point stored)
   
   properties (SetAccess = public, GetAccess = public, Hidden = true)
      ExtraQscore = [];    % The score matrix for extra quadratic constraints
      FeasibleFlag = false; % Whether the feasibility includes fitting error
   end
   
   properties
      Name                 % Descriptive name for the dataset
      DatasetUnits = [];   % B2BDC.B2Bdataset.DatasetUnitList object 
      Variables = [];      % B2BDC.B2Bvariables.VariableList object
      FeasiblePoint = [];  % Feasible point if the dataset is consistent
   end
   
   properties (Dependent)
      Length     % Number of dataset units in the dataset
   end
   
   properties (Dependent, Hidden = false)
      VarNames   % Cell array containing all variable names in the dataset
   end
   
   properties (SetAccess = private, GetAccess = public)
      ConsistencyMeasure = [];      % Lower and upper bounds to the consistency measure, [lowerBound upperBound]                
      ConsistencySensitivity = [];  % Structure array of sensitivity to consistency measure
      ModelDiscrepancy = []; % Information of model discrepancy correction
      ParameterDiscrepancy = []; % Information of parameter discrepancy correction
   end
   
   properties (SetAccess = private, GetAccess = public, Hidden = true)
      ModelDiscrepancyFlag = false; % Whether model discrepancy is considered in the analysis
      ParameterDiscrepancyFlag = false; % Whether parameter discrepancy is considered in the analysis
   end
   
   methods
      function obj = Dataset(dsName)
        %   OBJ = DATASET(DSNAME) returns a B2BDC.B2Bdataset.Dataset 
        %   object, where DSNAME is a descriptive name for the dataset.
        
         if nargin > 0
            if ischar(dsName)
               obj.Name = dsName;
            else
               error('Dataset name needs to be a character string')
            end
         end
         obj.Variables = B2BDC.B2Bvariables.VariableList();
      end
      
%       function b = loadobj(a)
%          % modified loading process of object
%          b = B2BDC.B2Bdataset.Dataset(a.Name);
%          b.DatasetUnits = a.DatasetUnits;
%          b.Variables = a.Variables;
%          b.ModelDiscrepancy = a.ModelDiscrepancy;
%          b.ParameterDiscrepancy = a.ParameterDiscrepancy;
%          b.FeasibleFlag = a.FeasibleFlag;
%          b.ExtraQscore = a.ExtraQscore;
%          b.ConsistencyMeasure = a.ConsistencyMeasure;
%          b.ConsistencySensitivity = a.ConsistencySensitivity;
%          b.ModelDiscrepancyFlag = a.ModelDiscrepancyFlag;
%          b.ParameterDiscrepancyFlag = a.ParameterDiscrepancyFlag;
%          b.FeasiblePoint = a.FeasiblePoint;
%       end
      
      function addDSunit(obj,dsUnitObj)
         %   ADDDSUNIT(OBJ, DSUNITOBJ) adds a dataset unit object DSUNITOBJ
         %   to the dataset OBJ. DSUNITOBJ can be a single DatasetUnit 
         %   object or an array of DatasetUnit objects to be added to the
         %   dataset.
         
         if isempty(obj.DatasetUnits)
            obj.DatasetUnits = B2BDC.B2Bdataset.DatasetUnitList();
         end
         if isa(dsUnitObj,'B2BDC.B2Bdataset.DatasetUnit')
            for i = 1:length(dsUnitObj)
               if ~isempty(obj)
                  allName = {obj.DatasetUnits.Values.Name};
                  [~,id] = intersect(allName, dsUnitObj(i).Name);
                  if ~isempty(id)
                     error(['The dataset unit ' allName{id} ' is already exsit in the dataset'])
                  end
                  sold = obj.DatasetUnits.Values(1).ScenarioParameter;
                  if ~isempty(sold)
                     snew = dsUnitObj(i).ScenarioParameter;
                     n1 = sold.Name;
                     n2 = snew.Name;
                     [~,~,id] = intersect(n1,n2,'stable');
                     if length(id) ~= length(n1)
                        error('Dataset units have different sets of scenario parameter!');
                     elseif ~isempty(id)
                        dsUnitObj(i).ScenarioParameter.Value = dsUnitObj(i).ScenarioParameter.Value(id);
                        dsUnitObj(i).ScenarioParameter.Name = dsUnitObj(i).ScenarioParameter.Name(id);
                     end
                  end
               end
               obj.DatasetUnits = add(obj.DatasetUnits,dsUnitObj(i));
               obj.Variables = obj.Variables.addList(dsUnitObj(i).VariableList);
            end
            obj.clearConsis;
            obj.clearModelDiscrepancy;
            obj.clearParameterDiscrepancy;
         else
            error(['A DatasetUnit object is required as an input. Use '...
                'generateDSunit to create a DatasetUnit object before adding'])
         end
      end
      
      function y = isempty(obj)
         %   Y = ISEMPTY(OBJ) returns a logical false if the Dataset OBJ
         %   contains DatasetUnits. A logical true is returned when OBJ
         %   contains no DatasetUnits.
         
         y = isempty(obj.DatasetUnits);
      end
      
      function y = get.Length(obj)
         y = obj.DatasetUnits.Length;
      end
      
      function y = get.VarNames(obj)
         datasetVar = obj.Variables.Values;
         n_variable = length(datasetVar);
         y = cell(n_variable,1);
         for i = 1:n_variable
            variable = datasetVar(i);
            y{i} = variable.Name;
         end
      end
      
      function y = length(obj)
         y = obj.Length;
      end
      
      function y = isConsistent(obj,opt)
         %   Y = ISCONSISTENT(OBJ) returns a logical true if the
         %   conistency measure of the dataset OBJ is positive.
         %
         %   Y = ISCONSISTENT(OBJ, OPT) the optional input of OPT can be
         %   used to modify the printed output or optimization criteria.
         %   OPT must be a B2BDC.Option object.
         
         if nargin > 1
            if ~isa(opt,'B2BDC.Option.Option')
               error(['Option input must be a B2BDC.Option object. ' ...
                   'Use generateOpt to create a B2BDC.Option object.'])
            end
         else
            opt = generateOpt;
         end
         flag1 = opt.ConsistencyMeasure;
         switch flag1
            case 'relative'
               if isempty(obj.ConsistencyMeasure)
                  if polytest(obj)
                     obj.polyConsisrel(opt);
                  elseif networktest(obj)
                     obj.networkConsisrel(opt);
                  else
                     obj.evalConsistencyrel(opt);
                  end
               end
            case 'absolute'
               if isempty(obj.ConsistencyMeasure)
                  if polytest(obj)
                     obj.polyConsisabs(opt);
                  elseif networktest(obj)
                     obj.networkConsisabs(opt);
                  else
                     obj.evalConsistencyabs(opt);
                  end
               end
            case 'DClab'
               if isempty(obj.ConsistencyMeasure)
                  obj.evalConsistencyDClab();
               end
         end
         if obj.ConsistencyMeasure(1) >= 0
            y = true;
         else
            y = false;
         end
         if opt.Display
            if y
               disp('The dataset is consistent')
            elseif obj.ConsistencyMeasure(2) <= 0
               disp('The dataset is inconsistent')
            else
               disp('The dataset consistency is undetermined')
            end
         end
      end
      
      function deletedUnits = deleteUnit(obj,id)
          %   DELETEDUNITS = DELETEUNIT(OBJ, ID) removes specified 
          %   DatasetUnits from the dataset OBJ and will return an
          %   n-by-1 array of deleted DatasetUnits in DELETEDUNITS. 
          %   ID will specify which DatasetUnits to delete by either a cell
          %   array of the DatasetUnit names or a vector of indicies of the
          %   n number of DatasetUnits to be removed. 
          
          deletedUnits = [];
          idx = [];
          if iscell(id)
             allName = {obj.DatasetUnits.Values.Name};
             [~,idx] = intersect(allName, id);
          elseif isa(id, 'B2BDC.B2Bdataset.DatasetUnit')
             allName = {obj.DatasetUnits.Values.Name};
             idNames = {id.Name};
             [~,idx] = intersect(allName, idNames);
          elseif isa(id, 'double') && all(id <= obj.DatasetUnits.Length) && all(id > 0)
             idx = id;
          end
          if isempty(idx)
             error(['The ID to be removed was not found ', ...
                'or exceeded the dataset length.'])
          end
          units = obj.DatasetUnits.Values;
          obj.clearDataset;
          if ~isempty(idx)
             deletedUnits = units(idx);
             units(idx) = [];
          end
          obj.addDSunit(units);
      end
      
      function changeBound(obj,newBD,idx)
          %   CHANGEBOUNDS(OBJ, NEWBD) returns a dataset OBJ with new
          %   bounds for all DatasetUnits. NEWDB is an nUnit-by-2 or by-3 
          %   matrix defining the LowerBound, UpperBound and ObservedValue 
          %   for all DatasetUnits. If ObservedValue is not provided the 
          %   mean of the provided LowerBound and UpperBound will be used. 
          %
          %   CHANGEBOUNDS(OBJ, NEWBD, IDX) returns a dataset OBJ with new
          %   bounds for DatasetUnits specified by IDX. NEWBD can be either
          %   n-by-2 or by-3 matrix defining the LowerBound, UpperBound
          %   and ObservedValue for the n DatasetUnits to be changed. 
          %   If ObservedValue is not provided the mean of the provided 
          %   LowerBound and UpperBound will be used. IDX is either a cell 
          %   array of n DatasetUnit names or a vector of length n of 
          %   indicies of DatasetsUnits to be calculated. 
          
          allUnitName = {obj.DatasetUnits.Values.Name};
          if nargin == 3
             if length(idx) ~= size(newBD,1)
                error('Inconsistent new bounds size with number of dataset units') 
             elseif iscell(idx)
                for i = 1:length(idx)
                   [~,id] = intersect(allUnitName, idx{i});
                   if isempty(id)
                      error(['The dataset unit ' idx{i} ' is not in the dataset'])
                   end
                   oldUnit = obj.DatasetUnits.Values(id);
                   dsName = obj.DatasetUnits.Values(id).Name;
                   newUnit = oldUnit.changeBounds(newBD(i,:));
                   obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
                end
                obj.clearConsis;
             elseif isnumeric(idx)
                for i = 1:length(idx)
                   oldUnit = obj.DatasetUnits.Values(idx(i));
                   dsName = obj.DatasetUnits.Values(idx(i)).Name;
                   newUnit = oldUnit.changeBounds(newBD(i,:));
                   obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
                end
                obj.clearConsis;
             else
                error('Wrong input index')
             end
          elseif nargin == 2 && size(newBD,1) == obj.Length
             for i = 1:obj.Length
                oldUnit = obj.DatasetUnits.Values(i);
                dsName = obj.DatasetUnits.Values(i).Name;
                newUnit = oldUnit.changeBounds(newBD(i,:));
                obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
             end
             obj.clearConsis;
          else
             error('Wrong number of input arguments')
          end   
      end
      
      function changeVarBound(obj,newBD,idx)
         v1 = obj.Variables;
         v2 = v1.changeBound(newBD,idx);
         obj.Variables = v2;
      end
      
      function dsUnit = findDSunit(obj,unitName)
         %   DSUNIT = FINDDSUNIT(OBJ, UNITNAME) looks for a DatasetUnit
         %   object from the dataset OBJ whose name matches the string
         %   UNITNAME. If no entry is found, DSUNIT will be returned
         %   empty.
         
         allName = {obj.DatasetUnits.Values.Name};
         id_unit = find(strcmp(allName,unitName),1);
         if ~isempty(id_unit)
            dsUnit = obj.DatasetUnits.Values(id_unit);
         else
            dsUnit = [];
         end
      end
      
      function dsNew = clone(obj)
         %   DSNEW = CLONE(OBJ) returns a deep copy (clone) of the dataset
         %   OBJ.
         
         dsName = obj.Name;
         dsNew = B2BDC.B2Bdataset.Dataset(dsName);
         units = obj.DatasetUnits.Values;
         n = obj.Length;
         for i = 1:n
            dsNew.addDSunit(units(i));
         end
         dsNew.FeasibleFlag = obj.FeasibleFlag;
         dsNew.FeasiblePoint = obj.FeasiblePoint;
         dsNew.ConsistencyMeasure = obj.ConsistencyMeasure;
         dsNew.ConsistencySensitivity = obj.ConsistencySensitivity;
      end
      
      function bounds = calBound(obj)
         %   BOUNDS = CALBOUND(OBJ) returns a matrix of size nUnit-by-2 of
         %   the eQOI bounds from the DatasetUnits in the dataset OBJ,
         %   where nUnit is the number of DatasetUnits.
         
         n = obj.Length;
         bounds = zeros(n,2);
         bounds(:,1) = [obj.DatasetUnits.Values.LowerBound]';
         bounds(:,2) = [obj.DatasetUnits.Values.UpperBound]';
      end
      
      function observes = calObserve(obj)
         %   OBSERVES = CALOBSERVE(OBJ) returns a vector containing all the
         %   ObservedValues from the DatasetUnits in the dataset OBJ.
         
         observes = [obj.DatasetUnits.Values.ObservedValue]';
      end
      
      function y = eval(obj, X, DSIdx)
         %   Y = EVAL(OBJ, X) evaluates all DatasetUnits surrogate models at
         %   X sampled points to produce a matrix Y of size nSample-by-nUnit
         %   of all the model's outputs. X is a matrix of size
         %   nSample-by-nVariable, where nSample is the number of samples
         %   to be evaluated and nVariable is the number of variables in the
         %   model.
         %
         %   Y = EVAL(OBJ, X, DSIDX) an input of DSIDX can be passed as
         %   either a a cell array of DatasetUnit names or a vector of
         %   indicies of the nUnit number of DatasetsUnits to be calculated.
         
         if nargin > 0
            if nargin < 3
               DSIdx = 1:obj.Length;
            end
            if obj.ModelDiscrepancyFlag
               nMD = obj.ModelDiscrepancy.Variables.Length;
            else
               nMD = 0;
            end
            if obj.ParameterDiscrepancyFlag
               nPD = obj.ParameterDiscrepancy.Variables.Length;
            else
               nPD = 0;
            end
            if size(X,2) ~= obj.Variables.Length && size(X,2) ~= obj.Variables.Length+nMD+nPD
               X = X';
            end
            if size(X,2) ~= obj.Variables.Length && size(X,2) ~= obj.Variables.Length+nMD+nPD
               error('The input number of variables does not match the number of variables in the dataset')
            elseif size(X,2) == obj.Variables.Length+nMD+nPD && nMD+nPD > 0
               y = obj.eval_with_discrepancy(X,DSIdx);
               return;
            end
            nSample = size(X,1);
            vObj = obj.Variables;
            if nargin == 2
               nUnit = obj.Length;
               y = zeros(nSample,nUnit);
               for i = 1:nUnit
                  tmpModel = obj.DatasetUnits.Values(i).SurrogateModel;
                  y(:,i) = tmpModel.eval(X, vObj);
               end
            elseif nargin == 3
               if isnumeric(DSIdx)
                  nUnit = length(DSIdx);
                  y = zeros(nSample, nUnit);
                  for i = 1:nUnit
                     tmpModel = obj.DatasetUnits.Values(DSIdx(i)).SurrogateModel;
                     y(:,i) = tmpModel.eval(X, vObj);
                  end
               elseif iscell(DSIdx)
                  nUnit = length(DSIdx);
                  allDS = {obj.DatasetUnits.Values.Name}';
                  y = zeros(nSample, nUnit);
                  for i = 1:nUnit
                     [~,dsID] = intersect(allDS, DSIdx{i});
                     if isempty(dsID)
                        error(['The dataset unit ' DSIdx{i} ' is not in the dataset'])
                     else
                        tmpModel = obj.DatasetUnits.Values(dsID).SurrogateModel;
                        y(:,i) = tmpModel.eval(X, vObj);
                     end
                  end
               else
                  error('Wrong input index type')
               end
            else
               error('Wrong number of input arguments')
            end
         end
      end
      
      function set.FeasiblePoint(obj,x0)
         if isempty(x0)
            obj.FeasiblePoint = [];
         else
            obj.FeasiblePoint = obj.findFeasiblePoint(x0);
         end
      end

      function clearConsis(obj)
         %   CLEARCONSIS(OBJ) removes all properties of consistency of
         %   the dataset OBJ
         
         obj.ConsistencyMeasure = [];
         obj.ConsistencySensitivity = [];
         obj.clearFeasiblePoint;
         obj.FeasibleFlag = false;
      end
      
      function clearModelDiscrepancy(obj)
         %  CLEARPARAMETERDISCREPANCY(OBJ) reset the model discrepancy to
         %  default (no correction)
         obj.ModelDiscrepancy = [];
         obj.ModelDiscrepancyFlag = false;
         obj.clearConsis;
      end
      
      function [sv,sname] = getScenarioParameter(obj)
         ss = [obj.DatasetUnits.Values.ScenarioParameter]';
         sname = ss(1).Name;
         sv = zeros(obj.Length,length(sname));
         for i = 1:obj.Length
            sv(i,:) = ss(i).Value;
         end         
      end
      
      function clearFeasiblePoint(obj)
         obj.FeasiblePoint = [];
         if obj.ModelDiscrepancyFlag
            obj.ModelDiscrepancy.FeasiblePoint = [];
         end
         if obj.ParameterDiscrepancyFlag
            obj.ParameterDiscrepancy.FeasiblePoint = [];
         end
      end
      
      function clearParameterDiscrepancy(obj)
         %  CLEARMODELDISCREPANCY(OBJ) reset the model discrepancy to
         %  default (no correction)
         obj.ParameterDiscrepancy = [];
         obj.ParameterDiscrepancyFlag = false;
         obj.clearConsis;
      end
      
      function makeSubset(obj,idx)
         % MAKESUBSET(OBJ,IDX) makes a new dataset object, containing only
         % dataset units specified by idx.
         
         if iscell(idx)
            name1 = {obj.DatasetUnits.Values.Name}';
            [~,~,idx] = intersect(idx,name1,'stable');
         end
         units = obj.DatasetUnits.Values;
         obj.clearDataset;
         if ~isempty(idx)
            remainUnits = units(idx);
         end
         obj.addDSunit(remainUnits);
      end
      
   end
   
   methods (Static, Hidden = true)
       xval = q2sample(Mgd,idx,H,Xvals,V);
   end
   
   methods (Static)
      [Qall,idx] = findConicHull(Q1,Q2)
      [Qall,idx] = approxConicHull(Q1,Q2,n)
   end
   
   methods (Hidden = true)
       function clearDataset(obj)
          %   CLEARDATASET(OBJ) returns a dataset OBJ whose dataset units
          %   and variable lists have been removed.
          
           obj.DatasetUnits = B2BDC.B2Bdataset.DatasetUnitList();
           obj.Variables = B2BDC.B2Bvariables.VariableList();
           obj.ConsistencyMeasure = [];
           obj.ConsistencySensitivity = [];
       end
       
       function y = polytest(obj)
           %   Y = POLYTEST(OBJ) returns a logical output true if all surrogate models are
           %   polynomial models and at least one of them is not quadratic model.
           
           y = false;
           dsUnits = obj.DatasetUnits.Values;
           for i = 1:length(dsUnits)
              testModel = dsUnits(i).SurrogateModel;
              if isa(testModel,'B2BDC.B2Bmodels.PolyModel')
                 y = true;
                 continue
              elseif ~isa(testModel,'B2BDC.B2Bmodels.QModel')
                 y = false;
                 break;
              end
           end
       end
       
       function xx0 = findFeasiblePoint(obj,x0)
          if size(x0,1) == 1
             x0 = x0';
          end
          nVar = obj.Variables.Length;
          if obj.ModelDiscrepancyFlag
             nMD = obj.ModelDiscrepancy.Variables.Length;
             xMD = obj.ModelDiscrepancy.FeasiblePoint;
             x0 = [x0; xMD];
          else
             nMD = 0;
          end
          if obj.ParameterDiscrepancyFlag
             nPD = obj.ParameterDiscrepancy.Variables.Length;
             xPD = obj.ParameterDiscrepancy.FeasiblePoint;
             x0 = [x0; xPD];
          else
             nPD = 0;
          end
          nx = length(x0);
          if nx == nVar+nMD+nPD
             if obj.isFeasiblePoint(x0')
                xx0 = x0(1:nVar);
                if nMD > 0
                   obj.ModelDiscrepancy.FeasiblePoint = x0(nVar+1:nVar+nMD);
                end
                if nPD > 0
                   obj.ParameterDiscrepancy.FeasiblePoint = x0(nVar+nMD+1:end);
                end
             else
                error('The input point is infeasible')
             end
          elseif nx == nVar
             [~,xnew] = obj.isFeasiblePoint(x0');
             if ~isempty(xnew)
                xx0 = x0;
                if nMD > 0
                   obj.ModelDiscrepancy.FeasiblePoint = xnew(nVar+1:nVar+nMD)';
                end
                if nPD > 0
                   obj.ParameterDiscrepancy.FeasiblePoint = xnew(nVar+nMD+1:end)';
                end
             else
                error('The input point is infeasible')
             end
          else
             error('The input point has a wrong dimension')
          end
       end
       
       
       evalConsistencyabs(obj,b2bopt)
       evalConsistencyDClab(obj)
       evalConsistencyrel(obj,b2bopt)
       [Qunits, Qx, Qextra, n_extra, extraIdx, L, idRQ, LBD] = getInequalQuad(obj,bds,frac)
       J = getJacobian(obj)
       directionSearch(obj,theta,x0,B2Bopt)
       y = eval_with_discrepancy(obj,X,DSIdx);
       [d,xopt] = calculateDistanceCurve(obj,x,opt,C)

       %CVX Functions
       [y,s] = cvxconsisabs(obj,yin,frac,abE)
       [y,s] = cvxconsisquadrel(obj,b2bopt,frac)
%        [yout,sensitivity] = obj.sedumiconsisquadrel_old(b2bopt, abE);
       [minout,minSensitivity,xs] = cvxminouterbound(obj,QOIobj,frac,abE,rflag)
       [maxout,maxSensitivity,xs] = cvxmaxouterbound(obj,QOIobj,frac,abE,rflag)
       [y,s] = cvxconsisrel(obj,yin,frac,abE)
       [y,s] = cvxconsisquadabs(obj,b2bopt,abE)
       
       % Sedumi Functions
       [y,s] = sedumiconsisabs(obj,opt,abE)
       [y,s] = sedumiconsisquadabs(obj,opt,abE)
       [y,s] = sedumiconsisrel(obj,yin,opt,abE)
       [y,s] = sedumiconsisquadrel(obj,opt,abE)
       [minout,minSensitivity] = sedumiminouterbound(obj,QOIobj,frac,abE,rflag)
       [maxout,maxSensitivity] = sedumimaxouterbound(obj,QOIobj,frac,abE,rflag)
       
       % nonlinear functions
       [Qmin, Qmax, ss, xOpt, alpha] = preQOIfmincon_minB(obj,q,disflag,rflag,b2bopt,alpha)
       [Qmin, Qmax, s, xOpt, abE] = preQOIopti(obj,QOIobj,disflag,rflag,b2bopt)
       [Qmin, Qmax, s, xOpt, abE] = preQOIfmincon(obj,q,disflag,rflag,b2bopt)
       [yin_result,s,xopt,abE,flag] = relCMfmincon(obj,disflag,b2bopt)
       [yin_result,s,xopt] = relCMfminconNN(obj,disflag,b2bopt)
       [yin_result,s,xopt,abE,flag] = relCMopti(obj,disflag,b2bopt)
       [yin_result,s,xopt,abE,flag] = absCMfmincon(obj,disflag,b2bopt)
       [yin_result,s,xopt] = absCMfminconNN(obj,disflag,b2bopt)
       [yin_result,s,xopt,abE,flag] = absCMopti(obj,disflag,b2bopt)
       
       % optimization
       
       
       % Matlab Functions
       y = addlistener(obj)
       y = delete(obj)
       y = findobj(obj)
       y = findprop(obj)
       y = notify(obj)
       y = le(obj)
       y = lt(obj)
       y = ne(obj)
       y = ge(obj)
       y = gt(obj)
       y = eq(obj) 
       
       
   end
   
end
