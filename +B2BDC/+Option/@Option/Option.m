classdef Option < handle
   % A B2BDC.Option.Option object that specifies the optimization crieteria.
   % The syntax of the input and the meaning of each property are in the
   % following:
   % opt = B2BDC.Option(PropName, PropValue ... )
   % ---------------------------------------------------------------------
   % ConsistencyMeasure:
   %   'relative' -  Calculate consistency measure in the way that
   %                  maximize              y
   %                 subject to  (1-y)*(LB(i)-OB(i)) < M(i)-OB(i)
   %                              M(i)-OB(i) < (1-y)*(UB(i)-OB(i))
   %                 OB(i) is the observation value of ith experiment
   %                 LB(i) is the lower bound of ith experiment
   %                 UB(i) is the upper bound of ith experiment
   %                 M(i) is the surrogate model of ith experiment
   %   'absolute' -  Calculate consistency measure in the way that 
   %                  maximize          y
   %                 subject to   LB(i) + y < M(i)
   %                              M(i) < UB(i) - y
   %   'DClab'    -  Calculate consistency measure from DClab package
   % ---------------------------------------------------------------------
   % ExtraLinFraction:
   %   1 ~ 100  -  Include the indicated number of fraction of extra
   %               linear pairs in outer bound calculation
   %      -1    -  Include extra linear pairs with influence factor greater
   %               than 5% of the most influential pair in outer bound
   %               calculation
   % ---------------------------------------------------------------------
   % TolConsis: The tolerance used in bi-sectional algorithm in consistency
   %            mearure outer bound calculation
   % ---------------------------------------------------------------------
   % Display:
   %   true  - Display the progress information in the optimization
   %           procedure
   %   false - Display no progress information in the optimization
   %           procedure
   % ---------------------------------------------------------------------
   % AddFitError:
   %   true  - The calculation in consistency measure and model prediction
   %           will add maximum absolute fitting error of each surrogate
   %           model to the corresponding QOI bound.
   %   false - The calculation in consistency measure and model prediction
   %           will not add maximum absolute fitting error of each surrogate
   %           model to the corresponding QOI bound.
   % ---------------------------------------------------------------------
   % SelfInconsisFlow:
   %    true - When doing the selfInconsisAnalysis, it will show up the
   %           sensitivity plot for every self-inconsistent dataset unit
   %           before showing the final stats
   %   false - Show final stats directly
   % ---------------------------------------------------------------------
   % Prediction:
   %   both - Both outer and inner bounds are returned
   %   inner - Only inner bounds are returned
   %   outer - Only outer bounds are returned
   % ---------------------------------------------------------------------
   % SampleOption:
   %   A B2BDC.Option.SampleOption object
   % ---------------------------------------------------------------------
   % POPOption:
   %   A B2BDC.Option.POPOption object
   % ---------------------------------------------------------------------
   % OptimOption:
   %   A B2BDC.Option.OptimOption object
   
   
   % Created: July 12, 2015    Wenyu Li
   %  Modified: Septempber 2, 2015   Wenyu Li (Data normalization option added)
   
   
   properties
      ConsistencyMeasure = [];
      ExtraLinFraction = [];
      TolConsis = [];
      Display = [];
      AddFitError = [];
      SelfInconsisFlow = [];
      Prediction = [];
      SampleOption = [];
      POPOption = [];
      OptimOption = [];
   end
   
   methods
      function  obj = Option(inputCell)
         % To generate a B2BDC.Option object
         p = {'ConsistencyMeasure','ExtraLinFraction','TolConsis','Display',...
            'AddFitError','SelfInconsisFlow','Prediction',...
            'SampleOption','POPOption','OptimOption'};
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
                        obj.ConsistencyMeasure = inputCell{2*i};
                     case 2
                        obj.ExtraLinFraction = inputCell{2*i};
                     case 3
                        obj.TolConsis = inputCell{2*i};
                     case 4
                        obj.Display = inputCell{2*i};
                     case 5
                        obj.AddFitError = inputCell{2*i};
                     case 6
                        obj.SelfInconsisFlow = inputCell{2*i};
                     case 7
                        obj.Prediction = inputCell{2*i};
                     case 8
                        obj.SampleOption = inputCell{2*i};
                     case 9
                        obj.POPOption = inputCell{2*i};
                     case 10
                        obj.OptimOption = inputCell{2*i};
                  end
               else
                  error('Invalid input property names')
               end
            end
         end
         if isempty(obj.ConsistencyMeasure)
            obj.ConsistencyMeasure = 'relative';
         end
         if isempty(obj.ExtraLinFraction)
               obj.ExtraLinFraction = 0;
         end
         if isempty(obj.TolConsis)
               obj.TolConsis = 1e-4;
         end
         if isempty(obj.Display)
               obj.Display = true;
         end
         if isempty(obj.AddFitError)
            obj.AddFitError = false;
         end
         if isempty(obj.SelfInconsisFlow)
            obj.SelfInconsisFlow = true;
         end
         if isempty(obj.Prediction)
            obj.Prediction = 'both';
         end
         if isempty(obj.SampleOption)
            obj.SampleOption = generateSampleOpt;
         end
         if isempty(obj.POPOption)
            obj.POPOption = generatePOPOpt;
         end
         if isempty(obj.OptimOption)
            obj.OptimOption = generateOptimOpt;
         end
      end
      
      function set.ConsistencyMeasure(obj,str1)
         c = {'relative','absolute','DClab'};
         if any(strcmp(c,str1))
            obj.ConsistencyMeasure = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.ExtraLinFraction(obj,num1)
         if (num1 >= 0 && num1 <= 100) || num1 == -1
            obj.ExtraLinFraction = num1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.Display(obj,logic1)
         if islogical(logic1)
            obj.Display = logic1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.AddFitError(obj,logic1)
         if islogical(logic1)
            obj.AddFitError = logic1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.SelfInconsisFlow(obj,logic1)
         if islogical(logic1)
            obj.SelfInconsisFlow = logic1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.Prediction(obj,str1)
         c = {'both','inner','outer'};
         if any(strcmp(c,str1))
            obj.Prediction = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.SampleOption(obj,opt1)
         if isa(opt1,'B2BDC.Option.SampleOption')
            obj.SampleOption = opt1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.POPOption(obj,opt1)
         if isa(opt1,'B2BDC.Option.POPOption')
            obj.POPOption = opt1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.OptimOption(obj,opt1)
         if isa(opt1,'B2BDC.Option.OptimOption')
            obj.OptimOption = opt1;
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

