function plotQOISensitivity(obj,QOIsens,LCnames)
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
outflag = isempty(QOIsens.Outer);
inflag = isempty(QOIsens.Inner);
if outflag
   s = QOIsens.Inner.max;
else
   s = QOIsens.Outer.max;
end
f = figure('NumberTitle','off','ToolBar','none',...
   'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
   'Name','QOI Min & Max Value Outer Bound Sensitivity Plot');
dsUnits = obj.DatasetUnits.Values;
n_unit = obj.Length;
n_variable = length(s.varl);
nL = length(s.linu);
if nargin < 3 && nL > 0
   LCnames = cell(nL,1);
   for i = 1:nL
      LCnames{i} = ['Linear constraint ' num2str(i)];
   end
else
   LCnames = {};
end
if length(LCnames) ~= nL
   error('Wrong input linear constraint name dimension')
end
varNames = obj.VarNames;
p1 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p2 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p3 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p4 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p5 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p6 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p7 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
p8 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
pall = [p1;p2;p3;p4;p5;p6;p7;p8];
for i = 1:8
   pall(i).Visible = 'off';
end
if outflag
   p4.Visible = 'on';
else
   p2.Visible = 'on';
end
bg = uibuttongroup('Visible','off','Units','normalized',...
   'Position',[0.02,0.1,0.1,0.2],'SelectionChangeFcn',@bgselection);
r1 = uicontrol('Parent',bg,'Style','radiobutton','String','All',...
   'Units','normalized','Position',[0.1,0.75,0.8,0.3],'FontSize',12);
r2 = uicontrol('Parent',bg,'Style','radiobutton','String','Ranked',...
   'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',12);
bg2 = uibuttongroup('Visible','off','Units','normalized',...
   'Position',[0.02,0.35,0.1,0.2],'SelectionChangeFcn',@bgselection2);
r3 = uicontrol('Parent',bg2,'Style','radiobutton','String','Inner',...
   'Units','normalized','Position',[0.1,0.75,0.8,0.3],'FontSize',12);
r4 = uicontrol('Parent',bg2,'Style','radiobutton','String','Outer',...
   'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',12);
bg3 = uibuttongroup('Visible','off','Units','normalized',...
   'Position',[0.02,0.6,0.1,0.2],'SelectionChangeFcn',@bgselection3);
r5 = uicontrol('Parent',bg3,'Style','radiobutton','String','Pred min',...
   'Units','normalized','Position',[0.1,0.75,0.8,0.3],'FontSize',12);
r6 = uicontrol('Parent',bg3,'Style','radiobutton','String','Pred max',...
   'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',12);

bg.Visible = 'on';
bg2.Visible = 'on';
bg3.Visible = 'on';
set(bg,'SelectedObject',r2);
set(bg2,'SelectedObject',r4);
set(bg3,'SelectedObject',r6);
if outflag
   set(r4,'Visible','off');
   set(bg2,'SelectedObject',r3);
elseif inflag
   set(r3,'Visible','off');
   set(bg2,'SelectedObject',r4);
end

if outflag
   %p3
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p3,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p3,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p4
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
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
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p4)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p4)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')

   %p7
   s = QOIsens.Inner.min;
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p7,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p7,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   %p8
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
   [sexp_sort,expidx] = sort(-sexp,'descend');
   [svar_sort,varidx] = sort(-svar,'descend');
   nexp_select = min(10,length(sexp_sort));
   nvar_select = min(10,length(svar_sort));
   sexp_high = -sexp_sort(1:nexp_select);
   svar_high = -svar_sort(1:nvar_select);
   expidx_high = expidx(1:nexp_select);
   varidx_high = varidx(1:nvar_select);
   Xexp = cell(nexp_select,1);
   Yexp = zeros(nexp_select,2);
   Xvar = cell(nvar_select,1);
   Yvar = zeros(nvar_select,2);
   for i = 1:nexp_select
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p8)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p8)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
elseif inflag
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p1,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p1
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p1,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p2
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
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
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p2)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p2)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')

   %p5
   s = QOIsens.Outer.min;
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p5,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p5,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   %p6
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
   [sexp_sort,expidx] = sort(-sexp,'descend');
   [svar_sort,varidx] = sort(-svar,'descend');
   nexp_select = min(10,length(sexp_sort));
   nvar_select = min(10,length(svar_sort));
   sexp_high = -sexp_sort(1:nexp_select);
   svar_high = -svar_sort(1:nvar_select);
   expidx_high = expidx(1:nexp_select);
   varidx_high = varidx(1:nvar_select);
   Xexp = cell(nexp_select,1);
   Yexp = zeros(nexp_select,2);
   Xvar = cell(nvar_select,1);
   Yvar = zeros(nvar_select,2);
   for i = 1:nexp_select
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p6)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p6)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
else
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p1,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p1
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p1,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p2
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
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
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p2)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p2)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
   %p3
   s = QOIsens.Inner.max;
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p3,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p3,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   %p4
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
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
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p4)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p4)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
   %p5
   s = QOIsens.Outer.min;
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p5,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p5,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   %p6
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
   [sexp_sort,expidx] = sort(-sexp,'descend');
   [svar_sort,varidx] = sort(-svar,'descend');
   nexp_select = min(10,length(sexp_sort));
   nvar_select = min(10,length(svar_sort));
   sexp_high = -sexp_sort(1:nexp_select);
   svar_high = -svar_sort(1:nvar_select);
   expidx_high = expidx(1:nexp_select);
   varidx_high = varidx(1:nvar_select);
   Xexp = cell(nexp_select,1);
   Yexp = zeros(nexp_select,2);
   Xvar = cell(nvar_select,1);
   Yvar = zeros(nvar_select,2);
   for i = 1:nexp_select
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p6)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p6)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
   %p7
   s = QOIsens.Inner.min;
   if length(s.expl) == 1
      s11 = [s.expl, s.expu;zeros(1,2)];
   else
      s11 = [s.expl, s.expu];
   end
   subplot(1,2,1,'Parent',p7,'Units','normalized',...
      'Position',[0.08 0.1 0.38 0.8]);
   h1 = barh(s11);
   h1(1).EdgeColor = 'blue';
   h1(1).FaceColor = 'blue';
   h1(2).EdgeColor = 'red';
   h1(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   title('Normalized sensitivity wrt dataset observation bounds')
   ylabel('Dataset observation index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   
   s12 = [s.varl, s.varu; s.linl, s.linu];
   subplot(1,2,2,'Parent',p7,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   title('Normalized sensitivity wrt linear constraint bounds')
   ylabel('Linear constraint index')
   xlabel('Sensitivity coefficient')
   ylim = length(s.varl)+nL+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   %p8
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu;s.linl;s.linu];
   [sexp_sort,expidx] = sort(-sexp,'descend');
   [svar_sort,varidx] = sort(-svar,'descend');
   nexp_select = min(10,length(sexp_sort));
   nvar_select = min(10,length(svar_sort));
   sexp_high = -sexp_sort(1:nexp_select);
   svar_high = -svar_sort(1:nvar_select);
   expidx_high = expidx(1:nexp_select);
   varidx_high = varidx(1:nvar_select);
   Xexp = cell(nexp_select,1);
   Yexp = zeros(nexp_select,2);
   Xvar = cell(nvar_select,1);
   Yvar = zeros(nvar_select,2);
   for i = 1:nexp_select
      if expidx_high(i) > n_unit
         id = expidx_high(i) - n_unit;
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,2) = sexp_high(i);
      else
         id = expidx_high(i);
         Xexp{i} = [dsUnits(id).Name ' (' int2str(id) ')'];
         Yexp(i,1) = sexp_high(i);
      end
   end
   for i = 1:nvar_select
      if varidx_high(i) <= n_variable
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      elseif varidx_high(i) <= 2*n_variable + nL
         id = varidx_high(i) - 2*n_variable;
         Xvar{i} = LCnames{id};
         Yvar(i,1) = svar_high(i);
      else
         id = varidx_high(i) - 2*n_variable - nL;
         Xvar{i} = LCnames{id};
         Yvar(i,2) = svar_high(i);
      end
   end
   Xexp = flipud(Xexp);
   Yexp = flipud(Yexp);
   Xvar = flipud(Xvar);
   Yvar = flipud(Yvar);
   
   subplot(2,1,1,'Parent',p8)
   h3 = barh(Yexp,min(1,0.1*nexp_select));
   h3(1).EdgeColor = 'blue';
   h3(1).FaceColor = 'blue';
   h3(2).EdgeColor = 'red';
   h3(2).FaceColor = 'red';
   legend('wrt DS observation LB','wrt DS observation UB','Location','best');
   set(gca,'yticklabel',Xexp,'YLim',[0 nexp_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt observation bounds')
   
   subplot(2,1,2,'Parent',p8)
   h4 = barh(Yvar,min(0.1*nvar_select,1));
   h4(1).EdgeColor = 'blue';
   h4(1).FaceColor = 'blue';
   h4(2).EdgeColor = 'red';
   h4(2).FaceColor = 'red';
   legend('wrt linear constraint LB','wrt linear constraint UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(2) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt linear constraint bounds')
end


   function bgselection(ss,dd)
      n1 = ss.SelectedObject.String;
      n2 = bg2.SelectedObject.String;
      n3 = bg3.SelectedObject.String;
      switch n1
         case 'All'
            if strcmp(n2,'Outer')
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p5.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p1.Visible = 'on';
               end
            else
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p7.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p3.Visible = 'on';
               end
            end
         case 'Ranked'
            if strcmp(n2,'Outer')
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p6.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p2.Visible = 'on';
               end
            else
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p8.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p4.Visible = 'on';
               end
            end
      end
   end

   function bgselection2(ss,dd)
      n1 = ss.SelectedObject.String;
      n2 = bg.SelectedObject.String;
      n3 = bg3.SelectedObject.String;
      switch n1
         case 'Inner'
            if strcmp(n2,'All')
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p7.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p3.Visible = 'on';
               end
            else
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p8.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p4.Visible = 'on';
               end
            end
         case 'Outer'
            if strcmp(n2,'All')
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p5.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p1.Visible = 'on';
               end
            else
               if strcmp(n3,'Pred min')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p6.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p2.Visible = 'on';
               end
            end
      end
   end

   function bgselection3(ss,dd)
      n1 = ss.SelectedObject.String;
      n2 = bg.SelectedObject.String;
      n3 = bg2.SelectedObject.String;
      switch n1
         case 'Pred min'
            if strcmp(n2,'All')
               if strcmp(n3,'Outer')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p5.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p7.Visible = 'on';
               end
            else
               if strcmp(n3,'Outer')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p6.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p8.Visible = 'on';
               end
            end
         case 'Pred max'
            if strcmp(n2,'All')
               if strcmp(n3,'Outer')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p1.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p3.Visible = 'on';
               end
            else
               if strcmp(n3,'Outer')
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p2.Visible = 'on';
               else
                  for j = 1:8
                     pall(j).Visible = 'off';
                  end
                  p4.Visible = 'on';
               end
            end
      end
   end
end
