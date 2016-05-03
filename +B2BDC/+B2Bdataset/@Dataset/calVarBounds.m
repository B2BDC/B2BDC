function xBound = calVarBounds(obj,index,opt)
%   XBOUND = CALVARBOUNDS(OBJ) returns a nVariable-by-2 matrix containing 
%   the outer bounds of all nVariables inside the dataset OBJ. The outer 
%   bounds are constrained by both the dataset units and variables inside 
%   the B2BDC. B2Bdataset.Dataset object OBJ. 
%
%   XBOUND = CALVARBOUNDS(OBJ, INDEX) returns a n-by-2 matrix containing 
%   the outer bounds of variables specified by INDEX. INDEX is either a 
%   cell array of variable names or a vector of indicies of the n number of 
%   variables to be calculated. 
%   
%   XBOUND = CALVARBOUNDS(OBJ, INDEX, OPT) the optional input of OPT can be 
%   used to modify the printed output or optimization criteria. OPT is a 
%   B2BDC.Option object.

%  Created: Nov 30, 2015     Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
end
nVar = obj.Variables.Length;
if index > nVar
    error('Index cannot exceed the number of variables in the dataset.');
end
vars = obj.Variables;
allName = obj.VarNames;
tep_target = zeros(nVar+1,1);
if isempty(index)
   id = 2:nVar+1;
elseif iscell(index)
   [~,~,id] = intersect(index, allName, 'stable');
   id = id+1;
elseif isnumeric(index)
   id = index+1;
else
   error('Wrong input index format')
end
xBound = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior Variable Bounds');
% end
for i = 1:length(id)
   tep_model = tep_target;
   tep_model(id(i)) = 1;
   tep_QOI = generateModel(tep_model,vars);
   QOIRange = obj.predictQOI(tep_QOI, opt);
   xBound(i,1) = QOIRange.min(1);
   xBound(i,2) = QOIRange.max(2);
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end