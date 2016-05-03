classdef PolyBasis
   % A PolyBasis object specifies the basis monomials used in a
   % polynomial function.
   
   % Created: Nov 12, 2015     Wenyu Li
   properties
      Value = []; % A nMonomial-by-nVar matrix defines the basis
   end
   
   properties (Dependent)
      Length  % Number of monomials in the basis
      Dimension  % Number of variables in the basis
   end
   
   methods
      function  obj = PolyBasis(pBasis)
         % Create a PolyBasis object.
         % The input argument is:
         %   - pBasis: A nMonomial-by-nVar matrix defines the basis. So
         %             each row r defines a monomial in the basis and the
         %             monomial is expressed as m = $\prod x_i^{r_i}$
         if nargin > 0
            if ismatrix(pBasis)
               obj.Value = pBasis;
            else
               error('The input basis should be a nBasis-by-nVar matrix')
            end
         end
      end
      
      function pBasis2 = plus(obj,pBasis1)
         % Returns a new PolyBasis that defined by the union of monomials from
         % the two component PolyBasis objects. The two PolyBasis should have the
         % same dimension.
         if isa(pBasis1,'B2BDC.B2Bvariables.PolyBasis')
            basis1 = obj.Value;
            basis2 = pBasis1.Value;
            if obj.Dimension ~= pBasis1.Dimension
               error('The added polynomial basis should have the same dimension')
            end
            newBasis = union(basis1,basis2,'rows');
            pBasis2 = B2BDC.B2Bvariables.PolyBasis(newBasis);
         else
            error('Wrong input class type');
         end
      end
      
      function  y = get.Length(obj)
         % Returns the length of the PolyBasis object.
         y = size(obj.Value,1);
      end
      
      function  y = get.Dimension(obj)
         % Returns the Dimension of the PolyBasis object.
         y = size(obj.Value,2);
      end
      
      function newBasis = expandDimension(obj,oldVar,newVar)
         % Returns a new PolyBasis object with expanded dimension. The old 
         % PolyBasis is based on the old VariableList and the new PolyBasis
         % is based on the new VariableList. The new VariableList should
         % have all variables included in the old VariableList.
         % The input arguments are:
         %    obj - The old PolyBasis
         %    oldVar - The VariableList defines the old PolyBasis
         %    newVar - The VariableList defines the new PolyBasis
         oldName = {oldVar.Values.Name};
         newName = {newVar.Values.Name};
         [~,~,id] = intersect(oldName, newName, 'stable');
         if length(id) ~= obj.Length
            error('The new basis should have all variables of the old basis')
         end
         pBasis = zeros(obj.Length, newVar.Length);
         pBasis(:,id) = obj.Value;
         newBasis = B2BDC.B2Bvariables.PolyBasis(pBasis);
      end
      
      function pBasis = mtimes(pBasis1,pBasis2)
         % Returns a PolyBasis object. The returned PolyBasis contains all
         % possible monomials made by production of two monomials, each of
         % them comes from one of the component PolyBasis respectively. The
         % two component PolyBasis should have the same dimension.
         % The input arguments are:
         %   pBasis1 - Component PolyBasis object
         %   pBasis2 - Component PolyBasis object
         if isa(pBasis1,'B2BDC.B2Bvariables.PolyBasis')
            if isa(pBasis2,'B2BDC.B2Bvariables.PolyBasis')
               basis1 = pBasis1.Value;
               basis2 = pBasis2.Value;
               if pBasis1.Dimension ~= pBasis2.Dimension
                  error('The added polynomial basis should have the same dimension')
               end
               newBasis = zeros(pBasis1.Length*pBasis2.Length,pBasis1.Dimension);
               for i = 1:pBasis2.Length
                  count = (i-1)*pBasis1.Length+1;
                  newBasis(count:count+pBasis1.Length-1,:) = basis1+repmat(basis2(i,:),pBasis1.Length,1);
               end
               newBasis = unique(newBasis,'rows');
               pBasis = B2BDC.B2Bvariables.PolyBasis(newBasis);
            elseif (round(pBasis2)-pBasis2) == 0
               basis1 = pBasis1.Value;
               newBasis = pBasis2*basis1;
               pBasis = B2BDC.B2Bvariables.PolyBasis(newBasis);
            end
         elseif isa(pBasis2,'B2BDC.B2Bvariables.PolyBasis')
            if isa(pBasis1,'B2BDC.B2Bvariables.PolyBasis')
               basis1 = pBasis1.Value;
               basis2 = pBasis2.Value;
               if pBasis1.Dimension ~= pBasis2.Dimension
                  error('The added polynomial basis should have the same dimension')
               end
               newBasis = zeros(pBasis1.Length*pBasis2.Length,pBasis1.Dimension);
               for i = 1:pBasis2.Length
                  count = (i-1)*pBasis1.Length+1;
                  newBasis(count:count+pBasis1.Length-1,:) = basis1+repmat(basis2(i,:),pBasis1.Length,1);
               end
               newBasis = unique(newBasis,'rows');
               pBasis = B2BDC.B2Bvariables.PolyBasis(newBasis);
            elseif (round(pBasis1)-pBasis1) == 0
               basis1 = pBasis2.Value;
               newBasis = pBasis1*basis1;
               pBasis = B2BDC.B2Bvariables.PolyBasis(newBasis);
            end
         else
            error('Wrong input class type');
         end
      end
      
      function y = isSubset(obj,poly2)
         % Returns a logical result. If the first input PolyBasis contains
         % all monomials in the second PolyBasis, the result is true.
         % Otherwise it is false.
         % The input arguments are:
         %   obj - A PolyBasis object.
         %   poly2 - A PolyBasis object.
         B1 = poly2.Value;
         B2 = obj.Value;
         y = true;
         for i = 1:poly2.Length
            id = intersect(B1(i,:),B2,'rows');
            if isempty(id)
               y = false;
               break
            end
         end
      end
      
      function newPoly = makeSubset(obj,idx)
         % Returns a new PolyBasis object that contains part of the
         % monomials in the old PolyBasis. The selection of the monomials
         % is specified by the input index.
         % The input arguments are:
         %   obj - The old PolyBasis object
         %   idx - An index array specifies the selected monomials.
         if ~isvector(idx)
            error('Wrong input index')
         end
         if max(idx) > obj.Length
            error('Index exceeds the length of the original polynomial basis')
         end
         B1 = obj.Value;
         B2 = B1(idx,:);
         newPoly = B2BDC.B2Bvariables.PolyBasis(B2);
      end
      
      function newPoly = deleteMonomial(obj,idx)
         % Returns a new PolyBasis with monomials specified by the input
         % index deleted from the old PolyBasis.
         % The input arguments are:
         %   obj - The old PolyBasis
         %   idx - An index array specifies the deleted monomials.
         if ~isvector(idx)
            error('Wrong input index')
         end
         if max(idx) > obj.Length
            error('Index exceeds the length of the original polynomial basis')
         end
         B1 = obj.Value;
         B1(idx,:) = [];
         newPoly = B2BDC.B2Bvariables.PolyBasis(B1);
      end
      
   end
   
   methods (Hidden = true)

   end
   
end

