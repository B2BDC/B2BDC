function plotConsistencySensitivity(obj)
%   PLOTCONSISTENCYSENSITIVITY(OBJ) displays a plot of the sensitivity to 
%   the consistency measure with respect to all the experimental and
%   variable bounds in a dataset OBJ. 

% Created: June 29, 2015    Wenyu Li
% Modified: December 8, 2015 Jim Oreluk: Changed Default plot to 'Ranked'.
% Will also automatically check for consistency if not already checked.

s = obj.ConsistencySensitivity;
if isempty(s)
    disp('Calculating consistency measure:')
    obj.isConsistent;
    obj.plotConsistencySensitivity;
else
   f = figure('NumberTitle','off','ToolBar','none',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name',['Consistency Sensitivity Plot of Dataset ' obj.Name]);
   dsUnits = obj.DatasetUnits.Values;
   n_unit = obj.Length;
   n_variable = length(s.varl);
   varNames = obj.VarNames;
   p1 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);
   p2 = uipanel('Parent',f,'Units','normalized','Position',[0.15,0,0.85,1]);

   p1.Visible = 'off';
   bg = uibuttongroup('Visible','off','Units','normalized',...
      'Position',[0.02,0.1,0.1,0.2],'SelectionChangeFcn',@bgselection);
   r1 = uicontrol('Parent',bg,'Style','radiobutton','String','All',...
      'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',12);
   r2 = uicontrol('Parent',bg,'Style','radiobutton','String','Ranked',...
      'Units','normalized','Position',[0.1,0.15,0.8,0.3],'FontSize',12);
   bg.Visible = 'on';
   set(bg,'SelectedObject',r2);  % Set ranked as default.

   
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
%    ylim = ceil(length(s.expl)/5)*5;
   ylim = length(s.expl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);

   
   s12 = [s.varl, s.varu];
   subplot(1,2,2,'Parent',p1,'Units','normalized',...
      'Position',[0.54 0.1 0.38 0.8]);
   h2 = barh(s12);
   h2(1).EdgeColor = 'blue';
   h2(1).FaceColor = 'blue';
   h2(2).EdgeColor = 'red';
   h2(2).FaceColor = 'red';
   legend('wrt variable LB','wrt variable UB','Location','best');
   title('Normalized sensitivity wrt variable bounds')
   ylabel('Variable index')
   xlabel('Sensitivity coefficient')
%    ylim = ceil(length(s.varl)/5)*5;
   ylim = length(s.varl)+2;
   set(gca,'YLim',[0,ylim])
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   
   sexp = [s.expl;s.expu];
   svar = [s.varl;s.varu];
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
      if varidx_high(i) > n_variable
         id = varidx_high(i) - n_variable;
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,2) = svar_high(i);
      else
         id = varidx_high(i);
         Xvar{i} = [varNames{id} ' (' int2str(id) ')'];
         Yvar(i,1) = svar_high(i);
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
   legend('wrt variable LB','wrt variable UB','Location','best');
   set(gca,'yticklabel',Xvar,'YLim',[0 nvar_select+1]);
   xlim = get(gca,'XLim');
   xlim(1) = 0;
   set(gca,'XLim',xlim);
   title('Normalized sensitivity wrt variable bounds')
end

   function bgselection(ss,dd)
      n1 = ss.SelectedObject.String;
      switch n1
         case 'All'
            p1.Visible = 'on';
            p2.Visible = 'off';
         case 'Ranked'
            p1.Visible = 'off';
            p2.Visible = 'on';
      end
   end
end