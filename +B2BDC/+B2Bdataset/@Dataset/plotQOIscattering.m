function [err,id] = plotQOIscattering(obj,x0,logFlag,opt)
% PLOTQOISCATTERING(OBJ,X0,LOGFLAG,OPT) plots how each QOI's prediction in the
% dataset compares to the experimental uncertainty.

%  Created: Nov 1, 2017     Wenyu Li

nQ = obj.Length;
if nargin < 3 || isempty(logFlag)
   logFlag = false(nQ,1);
end
if nargin<4
   opt = generateOpt;
end
y0 = obj.eval(x0);
y0(logFlag) = exp(y0(logFlag));
obs = obj.calObserve;
obs(logFlag) = exp(obs(logFlag));
bds = obj.calBound;
abE = zeros(nQ,1);
if opt.AddFitError
   for i = 1:nQ
      abE(i) = obj.DatasetUnits.Values(i).SurrogateModel.ErrorStats.absMax;
   end
end
bds = bds+[-abE abE];
bds(logFlag,:) = exp(bds(logFlag,:));
iN = bds(:,1) <=0;
bds(iN,1) = 0.01*obs(iN);
dub = bds(:,2)-obs;
dlb = bds(:,1)-obs;
dQ = y0'-obs;
err = dQ./obs;
eub = dub./obs;
elb = -dlb./obs;
if opt.Display
   figure('NumberTitle','off',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name',['QoI prediction scatter plot of ' obj.Name]);
   ee = errorbar(1:nQ,zeros(nQ,1),elb,eub,'o','MarkerSize',1e-5,'LineWidth',2);
   hold on
   plot(1:nQ,err,'o','MarkerSize',5.5,'MarkerFaceColor','r','MarkerEdgeColor','r');
   ee.Color = [169 169 169]/255;
   xlabel('QoI index')
   ylabel('Relative error')
   legend('Uncertainty bound','Model prediction')
   set(gca,'XLim',[0,nQ+1],'FontSize',17,'FontWeight','Bold','Position',[0.08,0.15,0.85,0.75],'LineWidth',2);
end
id = [find(err > eub); find(err < -elb)];


