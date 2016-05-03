function dsUnitList = checkSelfConsistency(obj,opt)
%   DSUNITLIST = CHECKSELFCONSISTENCY(OBJ) returns an array of
%   self-inconsistent DatasetUnits. A self-inconsistent DatasetUnit is one
%   that there is no point in the parameter space that can produce the
%   experimental observation within its uncertainty bounds. 
%
%   DSUNITLIST = CHECKSELFCONSISTENCY(OBJ, OPT) returns an array of
%   self-inconsistent DatasetUnits. OPT is a B2BDC.Option object which can
%   be used to alter the display output.
 
  % Created: Oct 20, 2015      Wenyu Li
  
  if nargin < 2
     opt = generateOpt;
     opt.Display = false;
  end
  dsUnitList = [];
  n = obj.Length;
  if opt.Display
     h = waitbar(0,'Start to check dataset unit self-consistency');
  end
  for i = 1:n
     dsTest = B2BDC.B2Bdataset.Dataset;
     dsUnit = obj.DatasetUnits.Values(i);
     dsTest.addDSunit(dsUnit);
     dsTest.isConsistent(opt);
     if dsTest.ConsistencyMeasure(2) <= 0
        dsUnitList = [dsUnitList; dsUnit];
     end
     if opt.Display
        waitbar(i/n,h,['Checking dataset unit self-consistency ' num2str(i) ' / ' num2str(n)]);
     end
  end
%   if opt.Display
%      waitbar(1,h,'Deleting self-inconsistent dataset units from dataset');
%   end
  if opt.Display
     close(h);
  end