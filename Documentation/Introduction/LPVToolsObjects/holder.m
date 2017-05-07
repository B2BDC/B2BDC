
%
%% Grid-based LPV Data Structures
% 
% A grid-based LPV model is a collection of linear models on a 
% gridded domain of parameter values.
% The core data structure for grid-based LPV models is the 
% |pss| (denoting parameter-varying state space model), 
% which stores the LPV system as a state space array
% defined on a finite, gridded domain.  As a simple example, consider an
% LPV system $S(\rho)$ that depends on a single scalar parameter $\rho$ in the
% domain $\rho \in [a,b]$.  The infrastructure requires the user to
% specify the domain with a finite grid, e.g. $N$ points in the interval
% $[a,b]$. The toolbox contains an |rgrid| data object to
% facilitate the creation and manipulation of multivariable
% parameter domains. The user must also specify the values of the state
% space system $S(\rho)$ at each point in this gridded domain.  The
% |pss| object stores the state-space array data using the
% standard MATLAB Control System Toolbox |ss| object.  Thus the
% |pss| can be viewed as the parameter-varying extension of the
% standard |ss| object. To summarize, the LPV system $S(\rho)$ is
% represented by a |pss| data object which stores the gridded
% domain and the array that defines the state-space data at each point
% in the domain. 
% 
% 
% The notions of parameter-varying matrices and parameter-varying
% frequency responses arise naturally to complement the |pss|
% objects.  LPV systems are time-varying and hence frequency responses
% can not be used to represent the system behavior as parameters vary.
% However frequency responses are useful to gain intuition about the
% system performance at fixed locations in the operating domain.
% LPVTools represents parameter varying matrices and frequency
% responses by |pmat| and |pfrd| data objects,
% respectively.  These two data objects are both stored as a data array
% defined on a gridded domain. 
% A |pmat| stores a |double| array, while a |pfrd| stores 
% an array of frequency responses (|frd| object in the Control System Toolbox). 
% Table 1 shows the relation between the core grid-based 
% LPVTools data objects (|pmat|,
% |pss|, |pfrd|) and existing MATLAB objects. The first
% row of the Table (``Nominal'') shows the basic MATLAB objects:
% matrices are |double| objects, state-space systems are
% |ss| objects, and frequency responses are |frd| objects.
% The third row of Table~\ref{fig:objects}
% (``Nominal Gridded LPV'') shows the corresponding core LPV objects.  
% The main point is that the
% (|pmat|, |pss|, |pfrd|) objects should be viewed
% as parameter-varying extensions of the standard MATLAB and
% Control Systems Toolbox objects (|double|, |ss|, |frd|).
% 
% The second row of the table (``Uncertain'') shows the equivalent
% objects used to represent uncertainty: uncertain matrices, state space
% systems, and frequency responses are represented by |umat|,
% |uss|, and |ufrd| objects, respectively.  These objects
% are part of the MATLAB Robust Control Toolbox.  The Robust Control
% Toolbox models the uncertainty as a perturbation $\Delta$ wrapped in
% feedback around a nominal part $M$, i.e. uncertainty is represented
% using a linear fractional transformation.  Real parametric, complex
% parametric, and unmodeled dynamic uncertainties can be modeled.  The
% fourth row of Table~\ref{fig:objects} (``Uncertain Gridded LPV'')
% shows the corresponding parameter-varying objects with uncertainty:
% uncertain parameter-varying matrices, state space systems, and
% frequency responses are represented by |upmat|, |upss|,
% and |upfrd| objects, respectively. These objects enable the integration
% of uncertainty into LPV models. 
%
%
%
%% LFT-based LPV Data Structures
% 
% A key component of the LFT-based LPVTools infrastructure 
% is the core LFT data structure object, referred to as a |tvreal|  
% (denoting a time-varying parameter). 
% The |tvreal| object is used to create a time-varying, real
% valued scalar object. The |tvreal| has a range, denoting
% the maximum and minimum value that the time-varying scalar
% can assume, and a rate-bound denoting the maximum
% and minimum rate of change of the time-varying scalar.
% The |tvreal| is used to model individual time-varying
% parameters, and construct parameter dependent LFT matrices
% and systems. 
% LPVTools represents LFT-based parameter varying matrices and state-space
% systems by |plftmat| and |plftss| data objects, respectively. 
% The |plftmat|, and |plftss| objects are constructed
% using |tvreal| elements, using a syntax that is a direct parallel 
% to the |ureal| syntax that is used to define 
% |umat| and |uss| objects in the Robust Control Toolbox.
% 
% 
% Both |plftmat| and |plftss| objects are stored as uncertain objects 
% (i.e. |umat| and |uss|, respectively) which in MATLAB consist of the 
% LFT ($M$, $\Delta_\rho$), with $\Delta_\rho$ being a block of |ureal| objects.
% Each object's constituent |tvreal| objects provide the additional information necessary to 
% handle the $\Delta_\rho$ block like a block of time-varying parameters.
% Uncertainty can be integrated into the |plftmat|, and |plftss| objects, 
% allowing these data objects to model system with, and without uncertainty.
% The |plftmat| and |plftss| objects should be viewed %xxx cut?
% as LFT-based parameter-varying extensions of the standard MATLAB,
% Control System Toolbox, and Robust Control Toolbox objects |double|, |ss|, |umat|, and |uss| as seen in 
% rows five ("Nominal LFT LPV") and six ("Uncertain LFT LPV") in Table 1}.