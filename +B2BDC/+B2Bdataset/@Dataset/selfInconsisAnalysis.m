function dsNew = selfInconsisAnalysis(obj,opt)
%   DSNEW = SELFINCONSISANALYSIS(OBJ) returns a dataset object with
%   all self-inconsistent dataset units removed. A statistical plot will be
%   returned showing frequency of bounds with high sensitivity to the
%   consistency measure as well as average sensitivity value for the 
%   dataset. 
%
%   DSNEW = SELFINCONSISANALYSIS(OBJ, OPT) returns a dataset object with
%   all self-inconsistent dataset units removed. OPT is an B2BDC.Option
%   object where OPT.SelfInconsisFlow can be set to a logical true, which 
%   will display the sensitivity to the consistency measure for all 
%   self-inconsistent dataset units before displaying the overall dataset 
%   statistics.
 
 %  Created: Oct 28, 2015     Wenyu Li
 
 dsNew = obj.clone;
 if nargin < 2
    opt = generateOpt;
    opt.Display = false;
    opt.ExtraLinFraction = -1;
 end
 self_incon_dsunits = dsNew.checkSelfConsistency(opt);
 allName = obj.VarNames;
 n1 = length(allName);
 nameList = cell(2*(n1+1),1);
 for i = 1:n1
    nameList{i} = [allName{i}, ' UB'];
    nameList{n1+i} = [allName{i}, ' LB'];
 end
 nameList{end-1} = 'Experiment UB';
 nameList{end} = 'Experiment LB';
 count = zeros(2*(n1+1),1);
 n = length(self_incon_dsunits);
 sVal = zeros(2*(n1+1),1);
 for i = 1:n
    d1 = self_incon_dsunits(i);
    dtest = generateDataset(d1.Name);
    dtest.addDSunit(d1);
    dtest.isConsistent(opt);
    s = dtest.ConsistencySensitivity;
    [s1,id] = sort([s.varu;-s.varl],'descend');
    varName = dtest.VarNames;
    nx = length(varName);
    id(s1 < 0.1) = [];
    id0 = id;
    for j = 1:length(id)
       if id(j) > nx
          id(j) = id(j)-nx;
       end
    end
    n2 = length(id);
    [~,id1,~] = intersect(allName,varName(id));
    for j = 1:n2
       if id0(j) > nx
          count(id1(j)+n1) = count(id1(j)+n1)+1;
          sVal(id1(j)+n1) = sVal(id1(j)+n1)+s1(j);
       else
          count(id1(j)) = count(id1(j))+1;
          sVal(id1(j)) = sVal(id1(j))+s1(j);
       end
    end
    if s.expu > 0.1
       count(end-1) = count(end-1)+1;
       sVal(end-1) = sVal(end-1)+s.expu;
    elseif -s.expl > 0.1
       count(end) = count(end)+1;
       sVal(end) = sVal(end)-s.expl;
    end
    if opt.SelfInconsisFlow
       dtest.plotConsistencySensitivity;
       waitfor(gcf);
    end
 end
 nameList(count==0) = [];
 sVal(count==0) = [];
 count(count==0) = [];
 sVal = sVal./count;
 n3 = length(nameList);
 f = figure('NumberTitle','off','ToolBar','none',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name','Self-inconsistent dataset units statistics');
 p1 = uipanel('Parent',f,'Units','normalized','Position',[0.18,0,0.82,1]);
 p2 = uipanel('Parent',f,'Units','normalized','Position',[0.18,0,0.82,1]);
 p2.Visible = 'off';
 bg = uibuttongroup('Visible','off','Units','normalized',...
    'Position',[0.02,0.1,0.15,0.2],'SelectionChangeFcn',@bgselection);
 r1 = uicontrol('Parent',bg,'Style','radiobutton','String','Number of apperance',...
    'Units','normalized','Position',[0.1,0.45,0.8,0.3],'FontSize',11);
 r2 = uicontrol('Parent',bg,'Style','radiobutton','String','Average sensitivity value',...
    'Units','normalized','Position',[0.1,0.15,0.8,0.3],'FontSize',11);
 bg.Visible = 'on';
 
 [c1,id] = sort(count);
 tag = nameList(id);
 subplot(1,1,1,'Parent',p1,'Units','normalized',...
      'Position',[0.12 0.1 0.72 0.8]);
 h1 = barh(c1);
 set(gca,'YTickLabel',tag)
 set(gca,'YTick',1:n3)
 set(gca,'YLim',[0,n3+0.5])
 set(gca,'Position',[0.2,0.1,0.75,0.8])
 title('Counts of rate constants or experiment bounds with sensitivity > 0.1')
 
 [c2,id] = sort(sVal);
 tag = nameList(id);
 subplot(1,1,1,'Parent',p2,'Units','normalized',...
      'Position',[0.12 0.1 0.72 0.8]);
 h2 = barh(c2);
 set(gca,'YTickLabel',tag)
 set(gca,'YTick',1:n3)
 set(gca,'YLim',[0,n3+0.5])
 set(gca,'Position',[0.2,0.1,0.75,0.8])
 title('Average sensitivity of rate constants or experiment bounds with sensitivity > 0.1')
 
 
   function bgselection(ss,dd)
      n1 = ss.SelectedObject.String;
      switch n1
         case 'Number of apperance'
            p1.Visible = 'on';
            p2.Visible = 'off';
         case 'Average sensitivity value'
            p1.Visible = 'off';
            p2.Visible = 'on';
      end
   end
end