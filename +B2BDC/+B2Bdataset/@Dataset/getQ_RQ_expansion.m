function [idall,Qall,Nall,Dall,APD,bPD] = getQ_RQ_expansion(obj)
% [IDALL,QALL,NALL,DALL] = GETQ_RQ_EXPANSION(OBJ) calculates the corresponding
% variable index, quadratic coefficient matrix, rational quadratic numerator 
% coefficient matrix and rational quadratic denumenator coefficient matrix for
% consistency measure or prediction calculation.

% Created: Sep 9, 2018     Wenyu Li


allVarnames = obj.VarNames;
n_variable = length(allVarnames);
units = obj.DatasetUnits.Values;
n_units = length(units);
idall = cell(n_units,1);
Nall = cell(n_units,1);
Dall = cell(n_units,1);
Qall = cell(n_units,1);
APD = [];
bPD = [];
if obj.ModelDiscrepancyFlag
   MDinfo = obj.ModelDiscrepancy;
   MDBasis = MDinfo.Basis;
   nmd = MDinfo.CorrectionDimension;
   GroupIndex = MDinfo.GroupIndex;
   nMD = sum(nmd);
   if obj.ParameterDiscrepancyFlag
      PDinfo = obj.ParameterDiscrepancy;
      PDBasis = PDinfo.Basis;
      Hv = PDinfo.VariableRange;
      npd = PDinfo.CorrectionDimension;
      nPD = sum(npd);
      for j = 1:n_units
         tmodel = units(j).SurrogateModel;
         [~,~,id1] = intersect(tmodel.VarNames,allVarnames,'stable');
         gid = GroupIndex(j);
         if gid > 0
            id2 = (sum(nmd(1:gid-1))+1:sum(nmd(1:gid)))';
         else
            id2 = [];
         end
         id3 = [];
         for k = 1:length(id1)
            if npd(id1(k)) ~= 0
               id3 = [id3;(sum(npd(1:id1(k)-1))+1:sum(npd(1:id1(k))))'];
            end
         end
         idall{j} = [id1;id2+n_variable;id3+n_variable+nMD];
         if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
            error('Discrepancy correction is currently not available for rational quadratic models');
         elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
            nc1 = length(id2);
            nc2 = length(id3);
            tCoef = tmodel.CoefMatrix;
            nm = size(tCoef,1);
            Coef = zeros(nm+nc1+nc2);
            Coef(1:nm,1:nm) = tCoef;
            if nc1 > 0
               Coef(1,nm+1:nm+nc1) = 0.5*MDBasis{j};
               Coef(nm+1:nm+nc1,1) = 0.5*MDBasis{j};
            end
            if nc2 > 0
               nn = npd(id1);
               tmpbasis = PDBasis{j};
               for k1 = 1:length(id1)
                  idv1 = id1(k1);
                  if ~isempty(tmpbasis{k1})
                     AA = zeros(2,n_variable+nMD+nPD);
                     aa = zeros(2,length(idall{j}));
                     bb = [Hv(idv1,2); -Hv(idv1,1)];
                     aa(:,k1) = [1 -1];
                     aa(:,nm+nc1+sum(nn(1:k1-1)):nm+nc1+sum(nn(1:k1))-1) = [tmpbasis{k1}; -tmpbasis{k1}];
                     AA(:,idall{j}) = aa;
                     APD = [APD; AA];
                     bPD = [bPD; bb];
                     tVec = tCoef(2:end,k1+1);
                     tmpQ = repmat(tVec,1,nn(k1)).*repmat(tmpbasis{k1},nm-1,1);
                     Coef(2:nm,nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1))) = tmpQ;
                     Coef(nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1)),2:nm) = tmpQ';
                     for k2 = 1:k1
                        if ~isempty(tmpbasis{k2})
%                            if k2 ~= k1
                              tmpQ = repmat(tCoef(k1+1,k2+1),nn(k1),nn(k2)).*repmat(tmpbasis{k1}',1,nn(k2)).*repmat(tmpbasis{k2},nn(k1),1);
                              Coef(nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1)),nm+nc1+sum(nn(1:k2-1))+1:nm+nc1+sum(nn(1:k2))) = tmpQ;
                              Coef(nm+nc1+sum(nn(1:k2-1))+1:nm+nc1+sum(nn(1:k2)),nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1))) = tmpQ';
%                            else
%                               tmpQ = repmat(tCoef(k1+1,k2+1),nn(k1),nn(k2)).*repmat(tmpbasis{k1}',1,nn(k2)).*repmat(tmpbasis{k2},nn(k1),1);
%                               Coef(nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1)),nm+nc1+sum(nn(1:k2-1))+1:nm+nc1+sum(nn(1:k2))) = tmpQ;
%                            end
                        end
                     end
                     tmpA = tCoef(1,k1+1)*tmpbasis{k1};
                     Coef(1,nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1))) = tmpA;
                     Coef(nm+nc1+sum(nn(1:k1-1))+1:nm+nc1+sum(nn(1:k1)),1) = tmpA;
                  end
               end
            end
            Qall{j} = Coef;
         end
      end
   else
      for j = 1:n_units
         tmodel = units(j).SurrogateModel;
         [~,~,id1] = intersect(tmodel.VarNames,allVarnames,'stable');
         gid = GroupIndex(j);
         if gid > 0
            id2 = (sum(nmd(1:gid-1))+1:sum(nmd(1:gid)))';
         else
            id2 = [];
         end
         idall{j} = [id1;id2+n_variable;];
         if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
            error('Model discrepancy correction is currently not available for rational quadratic models');
         elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
            nc1 = length(id2);
            tCoef = tmodel.CoefMatrix;
            nm = size(tCoef,1);
            Coef = zeros(nm+nc1);
            Coef(1:nm,1:nm) = tCoef;
            if nc1 > 0
               Coef(1,nm+1:nm+nc1) = 0.5*MDBasis{j};
               Coef(nm+1:nm+nc1,1) = 0.5*MDBasis{j};
            end
            Qall{j} = Coef;
         end
      end
   end
elseif obj.ParameterDiscrepancyFlag
   PDinfo = obj.ParameterDiscrepancy;
   PDBasis = PDinfo.Basis;
   Hv = PDinfo.VariableRange;
   npd = PDinfo.CorrectionDimension;
   nPD = sum(npd);
   for j = 1:n_units
      tmodel = units(j).SurrogateModel;
      [~,~,id1] = intersect(tmodel.VarNames,allVarnames,'stable');
      id3 = [];
      for k = 1:length(id1)
         if npd(id1(k)) ~= 0
            id3 = [id3;(sum(npd(1:id1(k)-1))+1:sum(npd(1:id1(k))))'];
         end
      end
      idall{j} = [id1;id3+n_variable];
      if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
         error('Parameter discrepancy correction is currently not available for rational quadratic models');
      elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
         nc2 = length(id3);
         tCoef = tmodel.CoefMatrix;
         nm = size(tCoef,1);
         Coef = zeros(nm+nc2);
         Coef(1:nm,1:nm) = tCoef;
         if nc2 > 0
            nn = npd(id1);
            tmpbasis = PDBasis{j};
            for k1 = 1:length(id1)
               idv1 = id1(k1);
               if ~isempty(tmpbasis{k1})
                  AA = zeros(2,n_variable+nPD);
                  aa = zeros(2,length(idall{j}));
                  bb = [Hv(idv1,2); -Hv(idv1,1)];
                  aa(:,k1) = [1 -1];
                  aa(:,nm+sum(nn(1:k1-1)):nm+sum(nn(1:k1))-1) = [tmpbasis{k1}; -tmpbasis{k1}];
                  AA(:,idall{j}) = aa;
                  APD = [APD; AA];
                  bPD = [bPD; bb];
                  tVec = tCoef(2:end,k1+1);
                  tmpQ = repmat(tVec,1,nn(k1)).*repmat(tmpbasis{k1},nm-1,1);
                  Coef(2:nm,nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1))) = tmpQ;
                  Coef(nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1)),2:nm) = tmpQ';
                  for k2 = 1:k1
                     if ~isempty(tmpbasis{k2})
%                         if k2 ~= k1
                           tmpQ = repmat(tCoef(k1+1,k2+1),nn(k1),nn(k2)).*repmat(tmpbasis{k1}',1,nn(k2)).*repmat(tmpbasis{k2},nn(k1),1);
                           Coef(nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1)),nm+sum(nn(1:k2-1))+1:nm+sum(nn(1:k2))) = tmpQ;
                           Coef(nm+sum(nn(1:k2-1))+1:nm+sum(nn(1:k2)),nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1))) = tmpQ';
%                         else
%                            tmpQ = repmat(tCoef(k1+1,k2+1),nn(k1),nn(k2)).*repmat(tmpbasis{k1}',1,nn(k2)).*repmat(tmpbasis{k2},nn(k1),1);
%                            Coef(nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1)),nm+sum(nn(1:k2-1))+1:nm+sum(nn(1:k2))) = tmpQ;
%                         end
                     end
                  end
                  tmpA = tCoef(1,k1+1)*tmpbasis{k1};
                  Coef(1,nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1))) = tmpA;
                  Coef(nm+sum(nn(1:k1-1))+1:nm+sum(nn(1:k1)),1) = tmpA;
               end
            end
         end
         Qall{j} = Coef;
      end
   end
else
   for j = 1:n_units
      tmodel = units(j).SurrogateModel;
      [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
      if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
         Nall{j} = tmodel.Numerator;
         Dall{j} = tmodel.Denominator;
      elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
         Qall{j} = tmodel.CoefMatrix;
      end
   end
end
infFlag = find(isinf(bPD));
bPD(infFlag) = [];
APD(infFlag,:) = [];