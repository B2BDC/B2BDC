function dsUnit = generateDSunit(dsName,modelObj,QoIBD,QoIval)
%generateDSunit   Creates a DatasetUnit object
%   DSUNIT = generateDSunit(NAME, MODELOBJ, EBOUNDS) creates a DatasetUnit
%   object DSUNIT, with internal name NAME.  MODELOBJ is a Model object
%   and BOUNDS is a 1-by-2 numeric array, containing asserted bounds on
%   the true value of the feature/QOI predicted by MODEL.  The dataset
%   unit represents constraints on the variables in MODELOBJ, using the
%   lower and upper bounds in EBOUNDS.   
%
%   DSUNIT = generateDSunit(NAME, MODELOBJ, EBOUNDS, QOIVAL) uses QOIVAL
%   to specify the observed value (optional).  The observed value should
%   fall between the bounds in EBOUNDS.
%
%   See also generateModelByFitting.

% Returns a B2BDC.B2Bdataset.DatasetUnit object. 
% The input arguments are:
%
%   dsName  - A string defines the name of the dataset unit
%  modelObj - A B2BDC.B2Bmodel.Model object, usually is the surrogate model
%             that maps the input-output relation of the original model
%    QoIBD  - A length 2 vector defines the lower and upper bounds of the
%             Quantity of interest(QoI)
%    QoIval - Optional, a scalar value defines the observed value of QoI.

% Created: Nov 1, 2015     Wenyu Li

if nargin == 4
   dsUnit = B2BDC.B2Bdataset.DatasetUnit(dsName,modelObj,QoIBD(1),QoIBD(2),QoIval);
elseif nargin == 3
   dsUnit = B2BDC.B2Bdataset.DatasetUnit(dsName,modelObj,QoIBD(1),QoIBD(2));
else
   error('Wrong number of inputs')
end