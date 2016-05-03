function plotVarRange(obj,index,opt)

 % Plot the conservative projection of the feasible set on a 2-D variable
 % space. It is conservative so the feasible set is a subset of the
 % projection. The 2-D variable space is specified by either an index
 % vector or a cell array of variable names.
 
 %  Created: Oct 24, 2015    Wenyu Li
 
 eps = 1e-3;
 if nargin < 2
    error('Not enough inputs');
 end
 if nargin < 3
    opt = generateOpt;
    opt.Display = false;
 end
 if ~obj.isConsistent(opt)
    disp('The dataset is inconsistent')
    return
 end
 varName = obj.VarNames;
 if length(index) ~= 2
    error('The input index should be of length 2');
 end
 if iscell(index)
    [~,~,id] = intersect(index, varName, 'stable');
 elseif isnumeric(index)
    id = index;
 else
    error('wrong input class');
 end
 varList = B2BDC.B2Bvariables.VariableList;
 for i = 1:2
    Modelvar = obj.Variables.Values(id(i));
    varList = varList.add(Modelvar);
 end
 z1 = zeros(36,2);   % outer bound
 z2 = zeros(36,2);   % inner bound
 h = waitbar(0,'progressing');
 for i = 1:36;
    theta = 5*i;
    Coef = zeros(3,3);
    Coef(1,2) = 0.5*cosd(theta);
    Coef(1,3) = 0.5*sind(theta);
    Coef = Coef+Coef';
%     Coef(1,1) = -x0(1)*cosd(theta)-x0(2)*sind(theta);
    QOI = B2BDC.B2Bmodels.QModel(Coef,varList);   
    Coef = zeros(3);
    Coef(1,2) = 0.5*sind(theta);
    Coef(1,3) = -0.5*cosd(theta);
    Coef = Coef+Coef';
%     Coef(1,1) = sind(theta)*x0(1)-cosd(theta)*x0(2);
    Q = B2BDC.B2Bmodels.QModel(Coef,varList);
    dsUnit = B2BDC.B2Bdataset.DatasetUnit('extra constraint',Q,-eps,eps,0);
    obj.addDSunit(dsUnit);
    if obj.isConsistent(opt)
       z1(i,1) = obj.sedumiminouterbound(QOI,opt.ExtraLinFraction);
       z1(i,2) = obj.sedumimaxouterbound(QOI,opt.ExtraLinFraction);
%        z1(i,1) = obj.cvxminouterbound(QOI,opt.ExtraLinFraction);
%        z1(i,2) = obj.cvxmaxouterbound(QOI,opt.ExtraLinFraction);
    else
       z1(i,1) = 0;
       z2(i,2) = 0;
    end
    obj.deleteUnit({'extra constraint'});
    waitbar(i/36,h);
 end
 x1 = zeros(73,1);
%  x2 = zeros(73,1);
 y1 = zeros(73,1);
%  y2 = zeros(73,1);
delete(h);
 for i = 1:36
    theta = i*5;
    x1(i) = cosd(theta)*z1(i,1);
    y1(i) = sind(theta)*z1(i,1);
%     x2(i) = cosd(theta)*z2(i,1);
%     y2(i) = sind(theta)*z2(i,1);
    x1(i+36) = cosd(theta)*z1(i,2);
    y1(i+36) = sind(theta)*z1(i,2);
%     x2(i+36) = cosd(theta)*z2(i,2);
%     y2(i+36) = sind(theta)*z2(i,2);
 end

%  id = (x1==0 & y1==0);
%  x1(id) = [];
%  y1(id) = [];
 
 x1(end) = x1(1);
%  x2(end) = x2(1);
 y1(end) = y1(1);
%  y2(end) = y2(1);

bds = obj.Variables.calBound;
rx = bds(id(1),:);
ry = bds(id(2),:);
x1(x1 < rx(1)) = rx(1);
x1(x1 > rx(2)) = rx(2);
y1(y1 < ry(1)) = ry(1);
y1(y1 > ry(2)) = ry(2);


%  plot(x2,y2,'bo','LineWidth',1.2)
% x3 = [-1,-1,1,1,-1];
% y3 = [-1,1,1,-1,-1];
% plot(x3,y3,'b-','LineWidth',1.2)

% hold on
f = figure('NumberTitle','off','ToolBar','none',...
   'Name','Correlated Variable Range Plot');
plot(x1,y1,'r-','LineWidth',1.2)
xlim(bds(id(1),:));
ylim(bds(id(2),:));
xlabel(varName{id(1)});
ylabel(varName{id(2)});

 