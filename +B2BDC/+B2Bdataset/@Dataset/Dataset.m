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
   
   properties
      Name                 % Descriptive name for the dataset
      DatasetUnits = [];   % B2BDC.B2Bdataset.DatasetUnitList object 
      Variables = [];      % B2BDC.B2Bvariables.VariableList object
      FeasiblePoint = [];  % Feasible point if the dataset is consistent
   end
   
   properties (Dependent)
      Length     % Number of dataset units in the dataset  
   end
   
   properties (Dependent, Hidden = true)
      VarNames   % Cell array containing all variable names in the dataset
   end
   
   properties (SetAccess = private, GetAccess = public)
      ConsistencyMeasure = [];      % Lower and upper bounds to the consistency measure, [lowerBound upperBound]                
   end
   
   properties (SetAccess = private, GetAccess = public, Hidden = true)
      ConsistencySensitivity = [];  % Structure array of sensitivity to consistency measure
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
               end
               obj.DatasetUnits = add(obj.DatasetUnits,dsUnitObj(i));
               obj.Variables = obj.Variables.addList(dsUnitObj(i).VariableList);
            end
            obj.clearConsis;
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
      
      function y = isConsistent(obj,opt)
         %   Y = ISCONSISTENT(OBJ) returns a logical true if the
         %   conistency measure of the dataset OBJ is positive.
         %
         %   Y = ISCONSISTENT(OBJ, OPT) the optional input of OPT can be
         %   used to modify the printed output or optimization criteria.
         %   OPT must be a B2BDC.Option object.
         
         if nargin > 1
            if ~isa(opt,'B2BDC.Option')
               error(['Option input must be a B2BDC.Option object. ' ...
                   'Use generateOpt to create a B2BDC.Option object.'])
            end
         else
            opt = B2BDC.Option();
         end
         flag1 = opt.ConsistencyMeasure;
         switch flag1
            case 'relative'
               if isempty(obj.ConsistencyMeasure)
                  if polytest(obj)
                     obj.polyConsisrel(opt)
                  else
                     obj.evalConsistencyrel(opt);
                  end
               end
            case 'absolute'
               if isempty(obj.ConsistencyMeasure)
                  obj.evalConsistencyabs(opt);
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
          if iscell(id)
             allName = {obj.DatasetUnits.Values.Name};
             [~,idx] = intersect(allName, id);
          else
             idx = id;
          end
          units = obj.DatasetUnits.Values;
          obj.clearDataset;
          if ~isempty(idx)
             deletedUnits = units(idx);
             units(idx) = [];
          end
          obj.addDSunit(units);
      end
      
      function changeBounds(obj,newBD,idx)
          %   CHANGEBOUNDS(OBJ, NEWBD) returns a dataset OBJ with new
          %   bounds for all DatasetUnits. NEWDB is an nUnit-by-2 or by-3 
          %   matrix defining the LowerBound, UpperBound and ObservedValue 
          %   for all DatasetUnits. If ObservedValue is not provided the 
          %   mean of the provided LowerBound and UpperBound will be used. 
          %
          %   CHANGEBOUNDS(OBJ, NEWBD, IDX) returns a dataset OBJ with new
          %   bounds for DatasetUnits specified by IDX. NEWDB can be either
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
            if size(X,2) ~= obj.Variables.Length;
               error('The input number of variables does not match the number of variables in the dataset')
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
               elseif iscell(DSIdx);
                  nUnit = length(DSIdx);
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
            if size(x0,1) == 1
               x0 = x0';
            end
            if obj.isFeasiblePoint(x0')
               obj.FeasiblePoint = x0;
            else
               error('The input point is infeasible')
            end
         end
      end
      
   end
   
   methods (Static, Hidden = true)
       xval = q2sample(Mgd,idx,H,Xvals,V);
   end
   
   methods (Hidden = true)
       
       function clearConsis(obj)
          %   CLEARCONSIS(OBJ) removes all properties of consistency of
          %   the dataset OBJ
          
           obj.ConsistencyMeasure = [];
           obj.ConsistencySensitivity = [];
           obj.FeasiblePoint = [];
       end
       
       function clearDataset(obj)
          %   CLEARDATASET(OBJ) returns a dataset OBJ whose dataset units
          %   and variable lists have been removed.
          
           obj.DatasetUnits = B2BDC.B2Bdataset.DatasetUnitList();
           obj.Variables = B2BDC.B2Bvariables.VariableList();
           obj.ConsistencyMeasure = [];
           obj.ConsistencySensitivity = [];
       end
       
       function y = polytest(obj)
           %   Y = POLYTEST(OBJ) returns a logical true if all surrogate
           %   models of the dataset OBJ are PolyModel. A logical false is
           %   returned if any surrogate model is not of the class
           %   PolyModel.
           
           y = true;
           for i = 1:obj.Length
               tep_model = obj.DatasetUnits.Values(i).SurrogateModel;
               if ~isa(tep_model,'B2BDC.B2Bmodels.PolyModel')
                   y = false;
                   break
               end
           end
       end
       
       evalConsistencyabs(obj,b2bopt)
       evalConsistencyDClab(obj)
       evalConsistencyrel(obj,b2bopt)
       [Qunits, Qx, Qextra, n_extra, extraIdx] = getInequalQuad(obj,bds,frac)
       J = getJacobian(obj)
       directionSearch(obj,theta,x0,B2Bopt)
       
       %CVX Functions
       [y,s] = cvxconsisabs(obj,yin,frac,tolerance)
       [y,s] = cvxconsisquadrel(obj,frac)
       [minout,minSensitivity] = cvxminouterbound(obj,QOIobj,frac)
       [maxout,maxSensitivity] = cvxmaxouterbound(obj,QOIobj,frac)
       [y,s] = cvxconsisrel(obj,yin,frac,tolerance)
       [y,s] = cvxconsisquadabs(obj,QOIobj,frac)
       
       % Sedumi Functions
       [y,s] = sedumiconsisabs(obj,opt,abE)
       [y,s] = sedumiconsisquadabs(obj,opt,abE)
       [y,s] = sedumiconsisrel(obj,yin,opt,abE)
       [y,s] = sedumiconsisquadrel(obj,opt,abE)
       [minout,minSensitivity] = sedumiminouterbound(obj,QOIobj,frac,abE)
       [maxout,maxSensitivity] = sedumimaxouterbound(obj,QOIobj,frac,abE)
       
       % Matlab Functions
       y = addlistener()
       y = delete()
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