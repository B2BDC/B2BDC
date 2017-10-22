function xVal = quadSample(Q,H,x0,V)
% XVAL = QUADSAMPLE(OBJ,Q,H,X0,V) returns a new feasible point. The
% feasibility is defined by quadratic inequalities saved in Q and box
% constraint in H. The direction is specified by V.

%  Created: Oct 18, 2017

if any(x0>H(:,2)) || any(x0<H(:,1))
   msgbox('The starting point is infeasible')
   xVal = [];
   return
else
   nQ = length(Q);
   for i1 = 1:nQ
      if [1;x0]'*Q{i1}*[1;x0] > 0
         msgbox('The starting point is infeasible')
         xVal = [];
         return
      end
   end
   [tp,tn] = findBoundary;
   tt = rand;
   tmove = tn+tt*(tp-tn);
   xVal = x0+tmove*V;
end

   function [tp,tn] = findBoundary
      hh = [H(:,2)-x0; x0-H(:,1)];
      vv = hh./[V;-V];
      tp = min(vv(vv>0));
      tn = max(vv(vv<0));
      for i = 1:nQ
         tmpQ = Q{i};
         c1 = tmpQ(1,1);
         v1 = 2*tmpQ(2:end,1);
         q1 = tmpQ(2:end,2:end);
         c = c1+x0'*v1+x0'*q1*x0;
         b = V'*v1+2*V'*q1*x0;
         a = V'*q1*V;
         delta = sqrt(b^2-4*a*c);
         r1 = (-b+delta)/2/a;
         r2 = (-b-delta)/2/a;
         if r1 > 0
            tp = min(tp,r1);
         end
         if r2 < 0
            tn = max(tn,r2);
         end
      end
   end
end