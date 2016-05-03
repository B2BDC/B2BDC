function xNew = expandBasis(x)
% Transform the data matrix x [x1,x2,...,xn] into an expanded basis data
% matrix [1,x1,x2,....,xn,x1^2,x1*x2,....,xn^2]. The input data matrix should
% be a tall matrix of size nSample-by-nVariable

% Created: July 21, 2015    Wenyu Li

[n_sample, n_variable] = size(x);
n_coef = 0.5*(n_variable+1)*(n_variable+2);
xNew = zeros(n_sample,n_coef);
xNew(:,1) = ones(n_sample,1);
xNew(:,2:n_variable+1) = x;
count = n_variable+2;
for i = 1:n_variable
   xNew(:,count:count+n_variable-i) = x(:,i:end).*repmat(x(:,i),1,n_variable+1-i);
   count = count+n_variable+1-i;
end