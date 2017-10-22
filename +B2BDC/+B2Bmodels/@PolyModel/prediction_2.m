function [min, max] = prediction(obj,N)
% Prediction of a polynomial model over its variable ranges. The input N
% defines the order of generalized Lagrangian multipliers, it should be a
% non-negative integer.

% Created: Nov 22, 2015     Wenyu Li

if nargin < 2
   N = 0;
end
vars = obj.Variables;
nVar = vars.Length;
lb = [vars.Values.LowerBound]';
ub = [vars.Values.UpperBound]';
pobj = obj.Basis;
cobj = obj.Coefficient;
p0 = pobj;
b0 = zeros(1,nVar);
for i = 1:nVar
   if lb(i)+ub(i) ~= 0
      b1 = repmat(b0,3,1);
      b1(2:3,i) = [1;2];
   else
      b1 = repmat(b0,2,1);
      b1(2,i) = 2;
   end
   p1 = B2BDC.B2Bvariables.PolyBasis(b1);
   b2 = repmat(b0,2*N+1,1);
   b2(:,i) = (0:2*N)';
   p2 = B2BDC.B2Bvariables.PolyBasis(b2);
   p0 = p0+p2*p1;
end
s = p0.sosPrune;
[s1,ids] = s.squareBasis;
if ~isSubset(s1,pobj)
   min = -inf;
   max = inf;
   return
end
[A,b,c,K,pars] = setSedumi();
[fopt,~,info] = sedumi(A,b,-c,K,pars);
min = fopt(1);
A(:,1) = -A(:,1);
b = -b;
[fopt,~,info] = sedumi(A,b,c,K,pars);
max = fopt(1);



   function [A,b,c,K,pars] = setSedumi()
      Li = (0:N)';
      p_tep = B2BDC.B2Bvariables.PolyBasis(Li);
      [sx, idx] = p_tep.squareBasis;
      nx = nVar*(N+1)^2+s.Length^2+1;
      K.f = 1;
      K.s = [repmat(N+1,1,nVar), s.Length];
      A = zeros(s1.Length, nx);
      b = zeros(s1.Length, 1);
      for i1 = 1:s1.Length
         b_tep = s1.Value(i1,:);
         id_tep = find(b_tep ~= 0);
         if isempty(id_tep)
            for i2 = 1:nVar
               c1 = -ub(i2)*lb(i2);
               if c1 ~= 0
                  [~,id] = intersect(sx.Value,0);
                  count = 1+(i2-1)*(N+1)^2;
                  id1 = idx{id}+count;
                  A(i1,id1) = c1*ones(1,length(id1));
               end
            end
            A(i1,1) = 1;
         elseif length(id_tep) == 1
            xd = b_tep(id_tep);
            count = 1+(id_tep-1)*(N+1)^2;
            [~,id] = intersect(sx.Value, xd-2);
            if ~isempty(id)
               id1 = idx{id}+count;
               A(i1,id1) = -1;
            end
            [~,id] = intersect(sx.Value, xd-1);
            if ~isempty(id)
               id1 = idx{id}+count;
               A(i1,id1) = ub(id_tep)+lb(id_tep);
            end
            [~,id] = intersect(sx.Value, xd);
            if ~isempty(id)
               id1 = idx{id}+count;
               A(i1,id1) = -ub(id_tep)*lb(id_tep);
            end
         end
         count = 1+nVar*(N+1)^2;
         id1 = count+ids{i1};
         A(i1,id1) = ones(1,length(id1));
         [~,id] = intersect(pobj.Value,b_tep,'rows');
         if ~isempty(id)
            b(i1) = cobj(id);
         end
         c = spalloc(nx,1,1);
         c(1) = 1;
         pars.fid = 0;
      end
      
   end




end
