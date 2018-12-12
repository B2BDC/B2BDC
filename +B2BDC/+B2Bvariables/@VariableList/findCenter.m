function x = findCenter(obj,x0,n)

% X=FINDCENTER(OBJ,X0,N) finds the center of the polytope defined by obj

if nargin < 2 || isempty(x0)
   x0 = obj.collectSamples(1);
   x0 = x0';
end
if nargin < 3
   n = 10;
end

fminopt1 = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','TolFun',...
   1e-8,'TolCon',1e-8);
   
% fminopt2 = optimoptions('fmincon','Display','none','GradObj','off',...
%    'GradConstr','on','Algorithm','interior-point','TolFun',...
%    1e-6,'TolCon',1e-6);

H = obj.calBound;
A = obj.ExtraLinConstraint.A;
if isempty(A)
   x = mean(H,2);
   return
end
LB = obj.ExtraLinConstraint.LB;
UB = obj.ExtraLinConstraint.UB;

[nC,nVar] = size(A);

dB = UB-LB;
A1 = [A eye(nC)];
A2 = [-A eye(nC)];
b = [UB;-LB];
xLB = [H(:,1); zeros(nC,1)];
xUB = [H(:,2); dB];

% 
% x1 = zeros(nVar,n);
% y1 = zeros(n,1);
% flag1 = zeros(n,1);
x2 = zeros(nVar,n);
y2 = zeros(n,1);
flag2 = zeros(n,1);

for i = 1:n
xS = [x0; dB.*rand(nC,1)];


% [xx,y1(i),flag1(i)] = fmincon(@funmax,xS,[A1;A2],b,[],[],...
%    xLB,xUB,[],fminopt1);
% x1(:,i) = xx(1:nVar);

[xx,y2(i),flag2(i)] = fmincon(@funmax_log,xS,[A1;A2],b,[],[],...
   xLB,xUB,[],fminopt1);
x2(:,i) = xx(1:nVar);

end
y2 = exp(-y2);
y2(flag2<0) = [];
x2(:,flag2<0) = [];
[~,id] = max(y2);
if ~isempty(id)
   x = x2(:,id)';
else
   x = [];
end



   function [y,gy] = funmax(x)
      y = -prod(x(nVar+1:nC+nVar));
      gy = -[zeros(nVar,1); y./x(nVar+1:nC+nVar)];
   end

   function [y,gy] = funmax_log(x)
      y = -sum(log(x(nVar+1:nC+nVar)));
      gy = -[zeros(nVar,1); 1./x(nVar+1:nC+nVar)];
   end


end