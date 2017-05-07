function [dsNew, changeID] = plotMinBoundChange(obj,opt)

 % Display a distribution of minimum experimental bounds change that can
 % bring the inconsistent dataset into consistent one. The change is
 % minimum such that the feasible set is almost(including some numerical effect)
 % a single point after the change.
 
 % Created: Oct 26, 2015     Wenyu Li

%  epsilon = 1e-5;

 if nargin < 2
    opt = generateOpt;
    opt.Display = false;
 end
 if obj.isConsistent(opt)
    disp('The dataset is already consistent')
    return
 end
 wy = obj.calObserve;
 [vec, flag] = obj.vectorConsistency(opt,[wy,wy]);
 dsNew = obj.clone;
 val = [obj.DatasetUnits.Values.ObservedValue]';
 dl = -vec.Relaxations.yLowerBound.*vec.Weights.yLowerBound;
 du = vec.Relaxations.yUpperBound.*vec.Weights.yUpperBound;
 changeID = union(find(-dl >= 1e-5), find(du >=1e-5));
 changeID = unique(changeID);
 bd1 = obj.calBound;
 bd2 = bd1;
 bd2(:,1) = bd2(:,1) + dl;
 bd2(:,2) = bd2(:,2) + du;
 dsNew.changeBound([bd2, val]);
 if opt.AddFitError
    dsNew.FeasibleFlag = true;
 end
 if flag == 1
    x0 = vec.FeasiblePoint;
    dsNew.FeasiblePoint = x0;
 end
 dl = dl./val;
 du = du./val;
 

 bd1(:,1) = (bd1(:,1)-val)./val;
 bd1(:,2) = (bd1(:,2)-val)./val;
 bd2(:,1) = (bd2(:,1)-val)./val;
 bd2(:,2) = (bd2(:,2)-val)./val;
 
 x = 1:obj.Length;
 y1 = zeros(obj.Length,1);
 bd3 = bd2;
 bd4 = bd2;
 bd3(:,1) = bd1(:,2);
 bd3(:,2) = bd2(:,2);
 bd4(:,1) = bd2(:,1);
 bd4(:,2) = bd1(:,1);
 y2 = mean(bd3,2);
 y3 = mean(bd4,2);
 e2 = 0.5*du;
 e3 = 0.5*dl;
 if opt.Display
    f = figure('NumberTitle','off','ToolBar','none',...
       'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
       'Name',['Bounds change of Dataset ' obj.Name]);
    a1 = axes('Parent',f);
    errorbar(x,y2,e2,'r.','LineWidth',1.2,'Parent',a1);
    hold on
    errorbar(x,y3,e3,'b.','LineWidth',1.2,'Parent',a1);
    errorbar(x,y1,bd1(:,1),bd1(:,2),'k.','LineWidth',1.5,'Parent',a1);
    
    title('The minimum bound change to make dataset consistent');
    legend('Relative change of UB','Relative change of LB','Original bounds')
    
    xlabel('Experiment index')
    ylabel('Relative uncertainty wrt observed experiment value')
    
    set(a1,'FontWeight','bold')
    set(a1,'FontSize',12)
    set(a1,'XLim',[0,obj.Length+1]);
 end