%% LPVTools Data Structures
% 
% 
% LPVTools is implemented using object-oriented programming. 
% The toolbox introduces several class-based data structures for modeling LPV systems. 
% These data structures extend the functionality associated with standard MATLAB data structures 
% from the Control Systems Toolbox and the Robust Control Toolbox into the LPV 
% framework. This is pictorially represented in Table 1. 
% 
% <<LPVToolsDataObjects.PNG>>
% 
% _Table 1: Relation between LPVTools and MATLAB objects._
%
% Table 1 shows the relation between the core 
% LPVTools data objects and existing MATLAB objects. 
% The first row of the Table (``Nominal'') shows the basic MATLAB objects:
% matrices are |double| objects, state-space systems are
% |ss| objects, and frequency responses are |frd| objects.
% The third row of Table~\ref{fig:objects}
% (``Nominal Gridded LPV'') shows the corresponding core grid-based LPV objects.  
% The core data structure for grid-based LPV models is the 
% |pss| (denoting parameter-varying state space model), 
% which stores the LPV system as a state space array
% defined on a finite, gridded domain.
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
% The (|pmat|, |pss|, |pfrd|) objects should be viewed
% as parameter-varying extensions of the standard MATLAB and
% Control Systems Toolbox objects (|double|, |ss|, |frd|).
%
% The second row of the table (``Uncertain'') shows the equivalent
% objects used to represent uncertainty: uncertain matrices, state space
% systems, and frequency responses are represented by |umat|,
% |uss|, and |ufrd| objects, respectively.  These objects
% are part of the MATLAB Robust Control Toolbox.
% The Robust Control
% Toolbox models the uncertainty as a perturbation $\Delta$ wrapped in
% feedback around a nominal part $M$, i.e. uncertainty is represented
% using a linear fractional transformation.  Real parametric, complex
% parametric, and unmodeled dynamic uncertainties can be modeled.  The
% fourth row of Table 1 (``Uncertain Gridded LPV'')
% shows the corresponding parameter-varying objects with uncertainty:
% uncertain parameter-varying matrices, state space systems, and
% frequency responses are represented by |upmat|, |upss|,
% and |upfrd| objects, respectively. These objects enable the integration
% of uncertainty into LPV models. 
% The (|upmat|, |upss|, |upfrd|) objects should be viewed
% as parameter-varying extensions of the uncertain Robust Control Toolbox 
% objects (|umat|, |uss|, |ufrd|).
%
% LPVTools represents LFT-based parameter varying matrices and state-space
% systems by |plftmat| and |plftss| data objects, respectively.
% Uncertainty can be integrated into the |plftmat|, and |plftss| objects, 
% allowing these data objects to model systems with, and without uncertainty.
% The |plftmat| and |plftss| objects should be viewed
% as LFT-based parameter-varying extensions of the standard MATLAB,
% Control System Toolbox, and Robust Control Toolbox objects |double|, |ss|, |umat|, and |uss|, as seen in 
% rows five ("Nominal LFT LPV") and six ("Uncertain LFT LPV") in Table 1}.


