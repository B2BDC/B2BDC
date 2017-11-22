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
   opt.Display = true;
   opt.ExtraLinFraction = 1;
end
pflag = opt.Prediction;
nVar = obj.Variables.Length;
vars = obj.Variables;
allName = obj.VarNames;
if isempty(index)
   id = 1:nVar;
elseif iscell(index)
   [~,~,id] = intersect(index, allName, 'stable');
elseif isnumeric(index)
   if any(index>nVar)
      error('Index cannot exceed the number of variables in the dataset.');
   end
   id = index;
else
   error('Wrong input index format')
end
selectName = allName(id);
xbd1 = zeros(length(id),2);
xbd2 = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior Variable Bounds');
% end
tep_model = [0;1];
for i = 1:length(id)
   tep_var = vars.makeSubset(id(i));
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
      
if opt.Display
   figure('NumberTitle','off',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name','Posterior parameter bounds');
   bds = vars.calBound;
   errorbar(id,zeros(length(id),1),abs(bds(id,1)),abs(bds(id,2)),'ok','MarkerSize',1e-5,'LineWidth',2);
   hold on
   switch pflag
      case 'inner'
         mbd = mean(xbd1,2);
         dx = xbd1(:,2)-mbd;
         errorbar(id-0.2,mbd,dx,'ob','MarkerSize',1e-5,'LineWidth',2);
         legend('Original variable bounds','B2B predicted inner bounds')
      case 'outer'
         mbd = mean(xbd1,2);
         dx = xbd1(:,2)-mbd;
         errorbar(id+0.2,mbd,dx,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original variable bounds','B2B predicted outer bounds')
      otherwise
         mbd1 = mean(xbd1,2);
         dx1 = xbd1(:,2)-mbd1;
         mbd2 = mean(xbd2,2);
         dx2 = xbd2(:,2)-mbd2;
         errorbar(id-0.2,mbd2,dx2,'ob','MarkerSize',1e-5,'LineWidth',2);
         errorbar(id+0.2,mbd1,dx1,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original variable bounds','B2B predicted inner bounds','B2B predicted outer bounds','LineWidth',2)
   end
   hold off
   xlabel('Variable index')
   ylabel('Uncertainty bounds')
   set(gca,'XLim',[0,max(id)+1],'FontSize',17,'FontWeight','Bold','Position',[0.08,0.15,0.85,0.75]);
end
