function quadCoef = vec2coef(a,n_variable)
% Transform the quadratic coefficient vector a into a matrix so the
% quadratic expression y = [1,x1,x2,...,xn,x1^2,x1*x2,...,xn^2] * a 
% into the expression y = [1,x1,x2,...,xn] * CoefMatrix * [1,x1,x2,...,xn]' 

% Created: July 21, 2015    Wenyu Li

n_Coef = 0.5*(n_variable+1)*(n_variable+2);
if ~isvector(a) || length(a) ~= n_Coef
   error('Wrong dimension of input coefficient or number of variables')
end
quadCoef = zeros(n_variable+1,n_variable+1);
count = 1;
for i = 1:n_variable+1
   quadCoef(i,i:end) = a(count:count+n_variable+1-i);
   count = count+n_variable+2-i;
end
quadCoef = 0.5*(quadCoef + quadCoef');