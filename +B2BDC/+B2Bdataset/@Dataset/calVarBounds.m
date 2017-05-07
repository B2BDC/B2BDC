function xPos = calVarBounds(obj,index,opt)
%   XPOS = CALVARBOUNDS(OBJ) returns a structure containing 
%   the posterior bounds of all variables inside the dataset OBJ. 
%
%   XPOS = CALVARBOUNDS(OBJ, INDEX) returns a structure containing 
%   the bounds of variables specified by INDEX. INDEX is either a 
%   cell array of variable names or a vector of indicies of the n number of 
%   variables to be calculated. 
%   
%   XBOUND = CALVARBOUNDS(OBJ, INDEX, OPT) the optional input of OPT can be 
%   used to modify the printed output or optimization criteria. OPT is a 
%   B2BDC.Option object.

%  Created: Nov 30, 2015     Wenyu Li
%  Modified: Nov 16, 2016    Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
   opt.ExtraLinFraction = -1;
end
pflag = opt.Prediction;
nVar = obj.Variables.Length;
if index > nVar
    error('Index cannot exceed the number of variables in the dataset.');
end
vars = obj.Variables;
allName = obj.VarNames;
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
selectName = allName(id-1);
xbd1 = zeros(length(id),2);
xbd2 = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior Variable Bounds');
% end
tep_model = [0;1];
for i = 1:length(id)
   tep_var = vars.makeSubset(id(i)-1);
   tep_QOI = generateModel(tep_model,tep_var);
   QOIRange = obj.predictQOI(tep_QOI, opt);
   if strcmp(pflag,'both')
      xbd1(i,1) = QOIRange.min(1);
      xbd1(i,2) = QOIRange.max(2);
      xbd2(i,1) = QOIRange.min(2);
      xbd2(i,2) = QOIRange.max(1);
   else
      xbd1(i,1) = QOIRange.min;
      xbd1(i,2) = QOIRange.max;
   end
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end
xPos.varName = selectName;
switch pflag
   case 'both'
      xPos.OuterBound = xbd1;
      xPos.InnerBound = xbd2;
   case 'inner'
      xPos.InnerBound = xbd1;
   case 'outer'
      xPos.OuterBound = xbd1;
end