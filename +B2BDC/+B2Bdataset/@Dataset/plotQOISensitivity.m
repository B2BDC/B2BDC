function plotQOISensitivity(obj,QOIsens)
% Plot sensitivity of QOI prediction with respect to experimental
% and variable bounds. The plot can display a general sensitivity
% distribution of all factors or a ranked plot including a few highest
% factors by selecting the corresponding button in the plot.
% The input argument is:
%   QOIsens - QOI sensitivity calculated from predictQOI function

% Created: August 4, 2015   Wenyu Li

if nargin < 2
   error('Not enough input')
end
smin = QOIsens.min;
smax = QOIsens.max;
f = figure('NumberTitle','off','ToolBar','none',...
   'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
   'Name','QOI Min & Max Value Outer Bound Sensitivity Plot');
dsUnits = obj.DatasetUnits.Values;
n_unit = obj.Length;
n_variable = length(smin.varl);
varNames = obj.VarNames;
p1 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p2 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
% p3 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p1.Visible = 'off';
% p3.Visible = 'off';
bg = uibuttongroup('Visible','off','Units','normalized',...
   'Position',[0.02,0.1,0.1,0.2],'SelectionChangeFcn',@bgselection);
r1 = uicontrol('Parent',bg,'Style','radiobutton','String','All',...
   'Units','normalized','Position',[0.1,0.75,0.8,0.3],'FontSize',12);
r2 = uicontrol('Parent',bg,'Style','radiobutton','String','Ranked',...
   'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',12);
% r3 = uicontrol('Parent',bg,'Style','radiobutton','String','Average',...
%    'Units','normalized','Position',[0.1,0.15,0.8,0.3],'FontSize',12);
bg.Visible = 'on';
set(bg,'SelectedObject',r2);

if length(smin.expl) == 1
   s11 = [smin.expl, smin.expu, smax.expl,  smax.expu;zeros(1,4)];
else
   s11 = [smin.expl, smin.expu, smax.expl,  smax.expu];
end
s1 = subplot(1,2,1,'Parent',p1,'Units','normalized',...
   'Position',[0.08 0.1 0.38 0.8]);
h1 = barh(s11,'Parent',s1);
h1(1).FaceColor = 'blue';
h1(1).EdgeColor = 'blue';
h1(2).FaceColor = 'red';
h1(2).EdgeColor = 'red';
h1(3).FaceColor = 'blue';
h1(3).EdgeColor = 'blue';
h1(4).FaceColor = 'red';
h1(4).EdgeColor = 'red';
l1 = legend('wrt exp LB','wrt exp UB','Location','northoutside',...
   'Orientation','horizontal');
ylabel('Dataset observation index')
xlabel('Sensitivity coefficient')
% ylim = ceil(length(smin.expl)/5)*5;
ylim = length(smin.expl)+2;
set(gca,'YLim',[0,ylim]);
text(0.05,0.95,'Predicted min','Units','normalized',...
   'HorizontalAlignment','left','FontSize',10,'FontWeight','bold')
text(0.95,0.95,'Predicted max','Units','normalized',...
   'HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
pp1 = l1.Position;
t1 = uicontrol('Style','text','String',...
   'Normalized sensitivity wrt dataset observation bounds',...
   'Units','normalized','FontWeight','bold',...
   'Position',[0.09,pp1(2)+0.03,0.36,0.03],...
   'FontSize',12,'Parent',p1);

s12 = [smin.varl, smin.varu, smax.varl,  smax.varu];
s2 = subplot(1,2,2,'Parent',p1,'Units','normalized',...
   'Position',[0.54 0.1 0.38 0.8]);
h2 = barh(s12,'Parent',s2);
h2(1).FaceColor = 'blue';
h2(1).EdgeColor = 'blue';
h2(2).FaceColor = 'red';
h2(2).EdgeColor = 'red';
h2(3).FaceColor = 'blue';
h2(3).EdgeColor = 'blue';
h2(4).FaceColor = 'red';
h2(4).EdgeColor = 'red';
l2 = legend('wrt var LB','wrt var UB','Location','northoutside',...
   'Orientation','horizontal');
ylabel('Variable index')
xlabel('Sensitivity coefficient')
% ylim = ceil(length(smin.varl)/5)*5;
ylim = length(smin.varl)+2;
set(gca,'YLim',[0,ylim])
text(0.05,0.95,'Predicted min','Units','normalized',...
   'HorizontalAlignment','left','FontSize',10,'FontWeight','bold')
text(0.95,0.95,'Predicted max','Units','normalized',...
   'HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
pp2 = l2.Position;
t2 = uicontrol('Style','text','String',...
   'Normalized sensitivity wrt variable bounds',...
   'Units','normalized','FontWeight','bold',...
   'Position',[0.55,pp2(2)+0.03,0.36,0.03],...
   'FontSize',12,'Parent',p1);


sexp = [-smin.expl; -smin.expu; smax.expl; smax.expu];
svar = [-smin.varl; -smin.varu; smax.varl; smax.varu];
[sexp_sort,expidx] = sort(sexp,'descend');
[svar_sort,varidx] = sort(svar,'descend');
nexp_select = min(10,length(sexp_sort));
nvar_select = min(10,length(svar_sort));
sexp_high = sexp_sort(1:nexp_select);
svar_high = svar_sort(1:nvar_select);
expidx_high = expidx(1:nexp_select);
varidx_high = varidx(1:nvar_select);
Xexp = cell(nexp_select,1);
Yexp = zeros(nexp_select,2);
Xvar = cell(nvar_select,1);
Yvar = zeros(nvar_select,2);
for i = 1:nexp_select
   if floor(0.5*(expidx_high(i)-1)/n_unit) == 0
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = -sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = -sexp_high(i);
      end
   else
      if expidx_high(i) > 3*n_unit
         id = expidx_high(i) - 3*n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i) - 2*n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
end
for i = 1:nvar_select
   if floor(0.5*(varidx_high(i)-1)/n_variable) == 0
      if varidx_high(i) > n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = -svar_high(i);
      else
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = -svar_high(i);
      end
   else
      if varidx_high(i) > 3*n_variable
         id = varidx_high(i) - 3*n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      end
   end
end
Xexp = flipud(Xexp);
Yexp = flipud(Yexp);
Xvar = flipud(Xvar);
Yvar = flipud(Yvar);

subplot(2,1,1,'Parent',p2)
h3 = barh(Yexp,min(0.8,0.1*nexp_select));
h3(1).EdgeColor = 'blue';
h3(1).FaceColor = 'blue';
h3(2).EdgeColor = 'red';
h3(2).FaceColor = 'red';
legend('wrt exp LB','wrt exp UB','Location','best');
set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1.5]);
title('Normalized sensitivity wrt dataset observation bounds')
xlabel('Sensitivity coefficient')
text(0.05,0.95,'Predicted min','Units','normalized',...
   'HorizontalAlignment','left','FontSize',10,'FontWeight','bold')
text(0.95,0.95,'Predicted max','Units','normalized',...
   'HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
xlim = get(gca,'XLim');
if xlim(2) > -xlim(1)
   xlim(1) = -xlim(2);
else
   xlim(2) = -xlim(1);
end
set(gca,'XLim',xlim);

subplot(2,1,2,'Parent',p2)
h4 = barh(Yvar,min(0.1*nvar_select,0.8));
h4(1).EdgeColor = 'blue';
h4(1).FaceColor = 'blue';
h4(2).EdgeColor = 'red';
h4(2).FaceColor = 'red';
legend('wrt var LB','wrt var UB','Location','best');
set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1.5]);
title('Normalized sensitivity wrt variable bounds')
xlabel('Sensitivity coefficient')
text(0.05,0.95,'Predicted min','Units','normalized',...
   'HorizontalAlignment','left','FontSize',10,'FontWeight','bold')
text(0.95,0.95,'Predicted max','Units','normalized',...
   'HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
xlim = get(gca,'XLim');
if xlim(2) > -xlim(1)
   xlim(1) = -xlim(2);
else
   xlim(2) = -xlim(1);
end
set(gca,'XLim',xlim);

% savg_exp = smax.expu-smax.expl+smin.expl-smin.expu;
% savg_var = smax.varu-smax.varl+smin.varl-smin.varu;
% savg = [savg_exp; savg_var];
% [sortAvg,id_avg] = sort(savg,'descend');
% navg_select = min(10,length(sortAvg));
% X = cell(navg_select,1);
% Y = zeros(navg_select,1);
% for i = 1:navg_select
%    if id_avg(i) > n_unit
%       id = id_avg(i) - n_unit;
%       X{i} = [varNames{id} ' (' int2str(id) ')'];
%    else
%       id = id_avg(i);
%       X{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
%    end
%    Y(i,1) = sortAvg(i);
% end
% X = flipud(X);
% Y = flipud(Y);
% subplot(1,1,1,'Parent',p3)
% h5 = barh(Y,min(0.1*nvar_select,0.8));
% h5(1).EdgeColor = 'blue';
% h5(1).FaceColor = 'blue';
% set(gca,'yticklabel',X,'YLim',[0 navg_select+0.5]);
% title('Ranked average sensitivity of prediction interval length wrt observation or variable bounds')
% xlabel('Sensitivity coefficient')
% xlim = get(gca,'XLim');
% xlim(1) = 0;
% set(gca,'XLim',xlim);







      function bgselection(ss,dd)
         n1 = ss.SelectedObject.String;
         switch n1
            case 'All'
               p1.Visible = 'on';
               p2.Visible = 'off';
               p3.Visible = 'off';
            case 'Ranked'
               p1.Visible = 'off';
               p2.Visible = 'on';
               p3.Visible = 'off';
            case 'Average'
               p1.Visible = 'off';
               p2.Visible = 'off';
               p3.Visible = 'on';
         end
      end
end
