function coefVector = coef2vec(obj)
% Transform the quadratic coefficient matrix into a vector coVec so the
% quadratic expression y = [1,x1,x2,...,xn] * CoefMatrix * [1,x1,x2,...,xn]' 
% into the expression y = [1,x1,x2,...,xn,x1^2,x1*x2,...,xn^2] * coefVector

% Created: July 21, 2015    Wenyu Li

quadCoef = 2*obj;
n_variable = size(quadCoef,1)-1;
n_coef = 0.5*(n_variable+1)*(n_variable+2);
coefVector = zeros(n_coef,1);
count = 1;
for i = 1:n_variable+1
   vec1 = quadCoef(i,i:end);
   vec1(1) = 0.5*vec1(1);
   coefVector(count:count+n_variable+1-i) = vec1;
   count = count+n_variable+2-i;
end