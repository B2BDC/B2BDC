function showSensitivity(X,Y,vars,ndis)

% SHOWSENSITIVITY(X,Y,VARS,NDIS) calculates and display the sensitivity of
% variables based on sample data X and Y.

[nS,nv] = size(X);
if vars.Length ~= nv
   error('Mismatched data and variable dimension');
end
if length(Y) ~= nS
   error('Mismatched data dimension');
end
xNew = [ones(size(X,1),1) X];
c = xNew \ Y;
xbd = vars.calBound;
dx = xbd(:,2) - xbd(:,1);
s = abs(c(2:end)).*dx;

figure('Name','Sensitivity Result','Units','normalized',...
   'Position',[0.2,0.1,0.5,0.7]);
if nargin < 4
   ndis = 15;
end
if nv < ndis
   ndis = nv;
end
[ss,id] = sort(s,'descend');
vName = {vars.Values.Name}';
s1 = ss(1:ndis);
n1 = vName(id(1:ndis));
for i = 1:ndis
   n1{i} = [n1{i} ' (' num2str(id(i)) ')'];
end
h = barh(flipud(s1));
h.BarWidth = 0.5;
aa = gca;
aa.YLim = [0,ndis+1];
aa.YTick = 1:ndis;
aa.YTickLabel = flipud(n1);
aa.FontSize = 13;
aa.FontWeight = 'bold';