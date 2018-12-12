function newOBJ = setScenarioParameter(oldOBJ,sv,sname)
% NEWOBJ = SETSCENARIOPARAMETER(OLDOBJ,SV,SNAME) resets the scenario
% parameter values to the input one.

% Created: Sep 9, 2018     Wenyu Li

if ~iscell(sname)
   sname = {sname};
end
if length(sv) ~= length(sname)
   error('Wrong input dimension');
else
   newOBJ = oldOBJ;
   newOBJ.ScenarioParameter.Value = sv;
   newOBJ.ScenarioParameter.Name = sname; 
end