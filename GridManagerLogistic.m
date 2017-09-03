classdef GridManagerLogistic < GridManagerLogExpAbstract
%This class is an implementation of the GridManager class for the Logistic
%model
%The meat of this class is in the GridManagerLogExpAbstract file
    
    properties (Constant)
        Name = 'Logistic';
        OverlappingGenerations = 1;
        ParamName1 = 'Birth Rate';
        ParamName2 = 'Death Rate';
        ParamBounds1 = [0 inf];
        ParamBounds2 = [0 inf];
        atCapacity = 0;
        plottingEnabled = 1;
    end
    
    methods (Access = public)
        
        function obj = GridManagerLogistic(dim, Ninit, mutation_manager, matrixOn, spatialOn, edgesOn, b, d)
            obj@GridManagerLogExpAbstract(dim, Ninit, mutation_manager, matrixOn, spatialOn, edgesOn, b, d);
        end

        function isSuccess = isSuccess(obj)
        	isSuccess = rand() > sum(obj.totalCount(:,obj.timestep))/obj.maxSize;
        end
    end
end