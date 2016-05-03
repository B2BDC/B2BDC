function [deleteUnit,dgap] = greedyDeletionFilter(dsOld,opt,factor)
 % Returns a dataset unit and the inner bound change after deleting 
 % this dataset unit which are selected to increase the inner bound of
 % consistency measure most from a set of dataset unit candidates. The
 % candidates come from the input factor that only dataset units with 
 % consistency measure sensitivity greater than or equal to factor will be
 % considered.
 % The input arguments are:
 %   dsOld - The dataset object
 %   opt - B2BDC.Option object
 %   factor - The threshold for dataset unit candidates (optional). If not
 %            specified, the default value is 0.05
 
 %  Created: Oct. 19th, 2015     Wenyu Li
 
 if nargin < 3
    factor = 0.05;
 end
 if nargin < 2
    opt = generateOpt;
    opt.Display = false;
 end
 if dsOld.isConsistent(opt);
    disp('The dataset is already consistent');
    deleteUnit = [];
    dgap = 0;
    return
 end
 d0 = dsOld.ConsistencyMeasure(1);
 s = dsOld.ConsistencySensitivity;
 n_Unit = dsOld.Length;
 allName = {dsOld.DatasetUnits.Values.Name};
 s_exp = [s.expu;-s.expl];
 [srank,idx] = sort(s_exp,'descend');
 id_opt = idx(srank>=factor);
 if isempty(id_opt)
    deleteUnit = {};
    dgap = 0;
 else
    n = length(id_opt);
    name_opt = cell(n,1);
    for i = 1:n
       if id_opt(i) > n_Unit
          id_opt(i) = id_opt(i)-n_Unit;
       end
       name_opt{i} = allName{id_opt(i)};
    end
    d_opt = zeros(n,1);
    optUnits = dsOld.DatasetUnits.Values(id_opt);
    for i = 1:n
       unitName = name_opt(i);
       dsUnit = optUnits(i);
       dsOld.deleteUnit(unitName);
       dsOld.isConsistent(opt);
       d_opt(i) = dsOld.ConsistencyMeasure(1);
       dsOld.addDSunit(dsUnit);
    end
    [dmax,i_unit] = max(d_opt);
    unitName = name_opt(i_unit);
    deleteUnit = dsOld.findDSunit(unitName);
    dgap = dmax-d0;
    dsOld.deleteUnit(unitName);
 end
 
 
       
       