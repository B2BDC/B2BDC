function ds = clearScenarioParameter(obj)
% ADDSCENARIOPARAMETER(OBJ,SVALUE,SNAME) adds the input scenario parameters 
% to the dataset units.

%  obj - A B2BDC.B2Bdataset.Dataset object
%  Svalue - A numerical array, each column correponds to one scenario parameter
%  Sname - A cell array, specifies which scenario parameters to be changed

dsUnits = obj.DatasetUnits.Values;
if ~isempty(obj.DatasetUnits.Values(1).ScenarioParameter.Value)
   for i = 1:obj.Length
      dsUnits(i).ScenarioParameter.Value = [];
      dsUnits(i).ScenarioParameter.Name = [];
   end
   ds = generateDataset(obj.Name);
   for i = 1:numel(dsUnits)
      ds.addDSunit(dsUnits(i));
   end
end