function varN = normalVar(vars,mx,sx)
 % Function used to normalize the VariableList object. It returns a
 % VariableList object with variable domain on a hypercube [-1, 1]^nVar
 % The input arguments are:
 %   mx - mean of variables
 %   sx - standard deviation of variables
 
 % Created: Oct 19, 2015      Wenyu Li

 allNames = {vars.Values.Name};
 n = vars.Length;
 H = [-ones(n,1), ones(n,1)];
 Val = [vars.Values.NominalValue];
 newVal = (Val-mx)./sx;
 varN = generateVar(allNames,H,newVal);