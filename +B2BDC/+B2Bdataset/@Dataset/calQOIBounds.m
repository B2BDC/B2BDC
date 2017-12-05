function qoiPos = calQOIBounds(obj, index, opt)
%   QOIPOS = CALQOIBOUNDS(OBJ) returns a structure containing 
%   the calculated bounds of all DatasetUnits inside the dataset OBJ,
%   where nDSUnit is the number of DatasetUnits inside the dataset.
%
%   QOIPOS = CALQOIBOUNDS(OBJ, INDEX) returns a structure containing 
%   the bounds of DatasetUnits specified by INDEX. INDEX is either a 
%   cell array of DatasetUnit names or a vector of indicies of the n number
%   of DatasetsUnits to be calculated. 
%   
%   QOIPOS = CALQOIBOUNDS(OBJ, INDEX, OPT) the OPT input can be used to 
%   modify the printed output or optimization criteria. OPT is a 
%   B2BDC.Option object.

%  Created: Nov 30, 2015     Wenyu Li
%  Modified: Nov 16, 2016     Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
   opt.ExtraLinFraction = 1;
end
pflag = opt.Prediction;
nUnit = obj.Length;
if index > nUnit
    error('Index cannot exceed the number of QOIs in the dataset.');
end
allName = {obj.DatasetUnits.Values.Name};
if isempty(index)
   id = 1:nUnit;
elseif iscell(index)
   [~,~,id] = intersect(index, allName, 'stable');
elseif isnumeric(index)
   id = index;
else
   error('Wrong input index format')
end
selectName = allName(id);
qbd1 = zeros(length(id),2);
qbd2 = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior QOI Bounds');
% end
for i = 1:length(id)
   tep_QOI = obj.DatasetUnits.Values(id(i)).SurrogateModel;
   QOIrange = obj.predictQOI(tep_QOI, opt);
   if strcmp(pflag,'both')
      qbd1(i,1) = QOIrange.min(1);
      qbd1(i,2) = QOIrange.max(2);
      qbd2(i,1) = QOIrange.min(2);
      qbd2(i,2) = QOIrange.max(1);
   else
      qbd1(i,1) = QOIrange.min;
      qbd1(i,2) = QOIrange.max;
   end
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end
qoiPos.varName = selectName;
switch pflag
   case 'both'
      qoiPos.OuterBound = qbd1;
      qoiPos.InnerBound = qbd2;
   case 'inner'
      qoiPos.InnerBound = qbd1;
   case 'outer'
      qoiPos.OuterBound = qbd1;
end

if opt.Display
   figure('NumberTitle','off',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name','Posterior QoI bounds');
   bds = obj.calBound;
   obs = obj.calObserve;
   dub = bds(:,2)-obs;
   dlb = obs-bds(:,1);
   errorbar(id,obs(id),dlb(id),dub(id),'ok','MarkerSize',1e-5,'LineWidth',2);
   hold on
   switch pflag
      case 'inner'
         mbd = mean(qbd1,2);
         dx = qbd1(:,2)-mbd;
         errorbar(id-0.2,mbd,dx,'ob','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted inner bounds')
      case 'outer'
         mbd = mean(qbd1,2);
         dx = qbd1(:,2)-mbd;
         errorbar(id+0.2,mbd,dx,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted outer bounds')
      otherwise
         mbd1 = mean(qbd1,2);
         dx1 = qbd1(:,2)-mbd1;
         mbd2 = mean(qbd2,2);
         dx2 = qbd2(:,2)-mbd2;
         errorbar(id-0.2,mbd2,dx2,'ob','MarkerSize',1e-5,'LineWidth',2);
         errorbar(id+0.2,mbd1,dx1,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted inner bounds','B2B predicted outer bounds')
   end
   hold off
   xlabel('QoI index')
   ylabel('Uncertainty bounds')
   set(gca,'XLim',[0,max(id)+1],'FontSize',17,'FontWeight','Bold','Position',[0.08,0.15,0.85,0.75]);
end