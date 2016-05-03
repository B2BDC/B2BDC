classdef DatasetUnitList < B2BDC.Util.IContainer
   % A B2BDC.B2Bdataset.DatasetUnitList is a container of 
   % B2BDC.B2Bdataset.DatasetUnit objects.
   
   % Created: June 15, 2015, myf + wl   
   
   methods
      
      function obj = DatasetUnitList()
         % Constructor of the B2BDC.B2Bdataset.DatasetUnitList object
         
         valClass = 'B2BDC.B2Bdataset.DatasetUnit';
         obj = obj@B2BDC.Util.IContainer(valClass);
      end
 
   end
   
end