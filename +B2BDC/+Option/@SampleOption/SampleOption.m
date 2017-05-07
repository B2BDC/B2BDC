classdef SampleOption < handle
   % A B2BDC.Option.SampleOption object that specifies the sampling crieteria.
   % The syntax of the input and the meaning of each property are in the
   % following:
   % opt = B2BDC.Option.SampleOption(PropName, PropValue ... )
   % ---------------------------------------------------------------------
   % SampleMethod:
   %   Sobol - Rejection sampling with sobolset on rotated box
   %   HR - Rejection sampling with hit and run in general polytope
   % ---------------------------------------------------------------------
   % UncertaintyEstimation:
   %   Outer - Outer bound of optimization
   %   Inner - Inner bound of optimization
   %   Sample - Sample based extreme estimation
   %   Sigma - Sample based 3 sigma estimation
   % ---------------------------------------------------------------------
   % BatchMaxSample:
   %   The maximum samples collected in one batch
   % ---------------------------------------------------------------------
   % StepInterval:
   %   The step interval in hit and run algorithm
   % ---------------------------------------------------------------------
   % PCTruncation:
   %   The threshold variance value for principal direction truncation
   % ---------------------------------------------------------------------
   % ExtraCut:
   %   The number of extra random direction cut used to generate containing
   %   polytope
   % ---------------------------------------------------------------------
   % UncertaintyTruncation:
   %   The relative truncation of the uncertainty bounds of polytope based
   %   on sampled uncertainty estimation
   % ---------------------------------------------------------------------
   % ParameterScaling:
   %   A logical value specifies whether scales parameter in sampling
   % ---------------------------------------------------------------------
   % PCinfo:
   %   A strcture specifying principal component information
   
   
   % Created: Oct 11, 2016    Wenyu Li
   %  Modified: Septempber 2, 2015   Wenyu Li (Data normalization option added)
   
   
   properties
      SampleMethod = [];
      UncertaintyEstimation = [];
      BatchMaxSample = [];
      StepInterval = [];
      PCTruncation = [];
      ExtraCut = [];
      UncertaintyTruncation = [];
      ParameterScaling = [];
      RejectionTol = [];
      PCinfo = [];
   end
   
   methods
      function  obj = SampleOption(inputCell)
         % To generate a B2BDC.Option object
         p = {'SampleMethod','UncertaintyEstimation','BatchMaxSample',...
            'StepInterval','PCTruncation','ExtraCut',...
            'UncertaintyTruncation','ParameterScaling',...
            'RejectionTol','PCinfo'};
         if nargin > 0 
            nin = length(inputCell);
         else
            nin = 0;
         end
         if mod(nin,2) ~= 0
            error('Wrong number of input argument')
         end
         nset = floor(0.5*nin);
         if nset > 0
            for i = 1:nset
               if any(strcmp(p,inputCell{2*i-1}))
                  id = find(strcmp(p,inputCell{2*i-1}) == 1);
                  switch id
                     case 1
                        obj.SampleMethod = inputCell{2*i};
                     case 2
                        obj.UncertaintyEstimation = inputCell{2*i};
                     case 3
                        obj.BatchMaxSample = inputCell{2*i};
                     case 4
                        obj.StepInterval = inputCell{2*i};
                     case 5
                        obj.PCTruncation = inputCell{2*i};
                     case 6
                        obj.ExtraCut = inputCell{2*i};
                     case 7
                        obj.UncertaintyTruncation = inputCell{2*i};
                     case 8
                        obj.ParameterScaling = inputCell{2*i};
                     case 9
                        obj.ParameterScaling = inputCell{2*i};
                     case 10
                        obj.ParameterScaling = inputCell{2*i};
                  end
               else
                  error('Invalid input property names')
               end
            end
         end
         if isempty(obj.SampleMethod)
            obj.SampleMethod = 'HR';
         end
         if isempty(obj.UncertaintyEstimation)
            obj.UncertaintyEstimation = 'Sample';
         end
         if isempty(obj.BatchMaxSample)
            obj.BatchMaxSample = 10^6;
         end
         if isempty(obj.StepInterval)
            obj.StepInterval = 10;
         end
         if isempty(obj.PCTruncation)
            obj.PCTruncation = 0;
         end
         if isempty(obj.ExtraCut)
            obj.ExtraCut = 0;
         end
         if isempty(obj.UncertaintyTruncation)
            obj.UncertaintyTruncation = 0;
         end
         if isempty(obj.ParameterScaling)
            obj.ParameterScaling = false;
         end
         if isempty(obj.RejectionTol)
            obj.RejectionTol = 1e-4;
         end
      end
      
      function set.SampleMethod(obj,str1)
         c = {'Sobol','HR'};
         if any(strcmp(c,str1))
            obj.SampleMethod = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.UncertaintyEstimation(obj,str1)
         c = {'Outer','Inner','Sample','Truncation'};
         if any(strcmp(c,str1))
            obj.UncertaintyEstimation = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.BatchMaxSample(obj,num1)
         c = round(num1);
         if c > 0
            obj.BatchMaxSample = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.StepInterval(obj,num1)
         c = round(num1);
         if c > 0
            obj.StepInterval = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.PCTruncation(obj,num1)
         if num1 >=0 && num1 < 1
            obj.PCTruncation = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.ExtraCut(obj,num1)
         c = round(num1);
         if c >= 0
            obj.ExtraCut = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.UncertaintyTruncation(obj,c)
         if c >= 0 && c < 1
            obj.UncertaintyTruncation = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.ParameterScaling(obj,c)
         if islogical(c)
            obj.ParameterScaling = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.RejectionTol(obj,c)
         if c > 0 && c < 1
            obj.RejectionTol = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.PCinfo(obj,s)
         if isstruct(s)
            names = fieldnames(s);
            if numel(names) == 3
               if strcmp(names{1},'direction') && strcmp(names{2},'variance') && strcmp(names{3},'mean')
                  obj.PCinfo = s;
               else
                  error('Invalid input property value')
               end
            else
               error('Invalid input property value')
            end
         else
            error('Invalid input property value')
         end
      end
   end
   
   methods (Hidden = true)
      % Matlab additional Functions
       y = addlistener() %matlab add listener events
       y = delete() % matlab delete handle
       y = findobj() % matlab find objects with a specified property value
       y = findprop() % matlab find matlab property of handle obj
       y = notify() 
       y = le() % matlab less than or equal to
       y = lt() % matlab less than
       y = ne() % matlab not equal to
       y = ge() % matlab greater than or equal to
       y = gt() % matlab greater than
       y = eq() % matlab equal to
%        y = isvalid() % matlab method for timer obj -- breaks when
       %hidden...
   end
   
end