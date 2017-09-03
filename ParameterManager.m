classdef ParameterManager < handle
%This powerhouse of a class stores the parameters for each of the different
%models and updates the display based on these parameters. It serves as an
%interface between the backend and frontend. 

    properties (Constant)
        maxTypes = 500; %The maximum possible number of different types 
        maxPopSize = 25000; %The maximum possible total population size (petri dish off)
        maxPlottingPopSize = 2500; %The maximum possible total population size (petri dish on)
        maxNumLoci = 15; %Maximum possible number of loci
    end

    properties (SetAccess = private)
    	matrixOn; %Determines whether the matrix model or vector model is used
        initialCounts; %Count of each genotype at inception
        handles; %Struct that stores pointers to the graphics objects
        classConstants; %Struct that stores the Constant properties of the currently used GridManager classes
        popSize; %The current total population size
        modelParameters; %The current parameters (Ninit, birth_rate, death_rate, fitness) of each type in each model
        numTypes; %The number of different types in the current model. NOT UPDATED BASED ON numLoci
        currentType; %The type who's parameters are currently displayed
        
        %Mutation 
        mutating; %Whether or not mutation is currently enabled.
        mutationMatrix; %The matrix that stores the current mutation transition probabilities
        recombining; %Whether or not recombination is enabled 
        recombinationNumber; %The recombination probability.
        
        %Multiple Loci
        s; %A parameter for determining the type parameters in the numLoci > 1 case
        numLoci; %The number of loci in each genotype
        e; %A parameter for determining the type parameters in the numLoci > 1 case
        
        %Non-Mutation Basic Parameters
        spatialOn; %Determines whether the spatial structure is considered.
        edgesOn; %Determines whether the edges are considered when searching for the closest square of a type. 
        currentModel; %The number of the current model

    end
    
    methods (Access = public)
        
        %Constructor method. This method initializes all of the instance
        %variables, particularly the seperate structs that store the
        %parameters for each model.
        function obj = ParameterManager(handles, classConstants)
            obj.classConstants = classConstants;
            obj.currentModel = 1;
            obj.currentType = 1;
            obj.setMatrixOn(handles.matrixOn_button.Value);
            obj.spatialOn = 1;
            obj.edgesOn = 1;
            %numerical parameters
            obj.popSize = 2500;
            obj.numTypes = 2;
            
            obj.handles = handles;
            obj.modelParameters = struct(...
                'Ninit_default', repmat({[]},1,length(obj.classConstants)),...
                'Param1_default', repmat({[]},1,length(obj.classConstants)),...
                'Param2_default', repmat({[]},1,length(obj.classConstants)),...
                'Ninit', repmat({[]},1,length(obj.classConstants)));
            for i = 1:length(obj.classConstants)
                if obj.classConstants(i).atCapacity
                    obj.modelParameters(i).Ninit_default = 1250;
                else
                    obj.modelParameters(i).Ninit_default = 1;
                end
                obj.modelParameters(i).Param1_default = 1;
                obj.modelParameters(i).Param2_default = 0.01;
                obj.modelParameters(i).Ninit = repmat(obj.modelParameters(i).Ninit_default,1, 2);
                obj.modelParameters(i).Param1 = repmat(obj.modelParameters(i).Param1_default,1, 2);
                obj.modelParameters(i).Param2 = repmat(obj.modelParameters(i).Param2_default,1, 2);
            end
            %mutations
            obj.mutating = 0;
            obj.mutationMatrix = [0.99 0.01; 0.01 0.99];
            obj.initialCounts = [1 0];
            obj.numLoci = 1;
            obj.recombining = 0;
            obj.recombinationNumber = 1;

            %multiple loci params
            obj.s = -0.5;
            obj.e = 0;
        end
        
        
        function writeBoxes(obj)
            %Writes the parameter manager's values to the box values
            %First, change the current type and current model
            obj.updateButtonValues();
            %Then, write to the boxes
            obj.handles.recombination_box.String = obj.recombinationNumber;
            obj.handles.loci_box.String = obj.numLoci;
            obj.handles.population_box.String = obj.popSize;
            if obj.mutating && obj.numLoci > 1
                 obj.handles.param_1_box.String = obj.s;
            else
                 obj.handles.param_1_box.String = obj.modelParameters(obj.currentModel).Param1(obj.currentType);
            end
            if obj.mutating && obj.numLoci > 1
                 obj.handles.param_2_box.String = obj.e;
            else
                 obj.handles.param_2_box.String = obj.modelParameters(obj.currentModel).Param2(obj.currentType);
            end
            obj.handles.init_pop_box.String = obj.modelParameters(obj.currentModel).Ninit(obj.currentType);
            obj.handles.num_types_box.String = obj.numTypes; 
        end

        
        function updateButtonValues(obj)
            %Called to update the values that are stored as yes/no buttons
            obj.currentType = obj.handles.types_popup.Value;
            obj.mutating = obj.handles.genetics_button.Value;
            obj.recombining = obj.handles.recombination_check.Value;
            obj.spatialOn = ~obj.handles.spatial_structure_check.Value;
            obj.matrixOn = obj.handles.matrixOn_button.Value;
            obj.edgesOn = ~obj.handles.remove_edges_check.Value;
            %current model
            newModel =  1*obj.handles.model1_button.Value + ...
                        2*obj.handles.model2_button.Value + ...
                        3*obj.handles.model3_button.Value + ...
                        4*obj.handles.model4_button.Value;
            if newModel ~= obj.currentModel
                obj.currentModel = newModel;
                obj.initializeInitialCounts(obj.numLoci, obj.popSize);
            end
        end
        
        function message = updateBoxValues(obj)
            %Updates the mutation/recombination and basic parameters to the
            %input values
            %Recombination, numTypes, param1, param2, initPop, maxPop, numLoci
            message = '';
            tempRecombinationNumber = str2double(obj.handles.recombination_box.String);
            tempNumLoci = str2double(obj.handles.loci_box.String);
            tempPopulation = str2double(obj.handles.population_box.String);
            tempParam1 = str2double(obj.handles.param_1_box.String);
            tempParam2 = str2double(obj.handles.param_2_box.String);
            ninitTemp = str2double(obj.handles.init_pop_box.String);
            numTypesTemp = str2double(obj.handles.num_types_box.String);
            if length(tempRecombinationNumber) ~= 1 || length(tempNumLoci) ~= 1 ||...
                            length(tempPopulation) ~= 1 || length(tempParam1) ~= 1 ||...
                            length(tempParam2) ~= 1 || length(ninitTemp) ~= 1 ||...
                            length(numTypesTemp) ~= 1
                message = sprintf('ERROR: All inputs must be scalar\n');
            end
            %Recombination Number
            if ParameterManager.isNumberWithin(tempRecombinationNumber, 0, 1)
                obj.recombinationNumber = tempRecombinationNumber;
            else
                message = sprintf('ERROR: Recombination Parameter must be a positive integer between 0 and 1\n');
            end
            %Number Loci
            if ParameterManager.isPositiveInteger(tempNumLoci)
                obj.updateNumLoci(tempNumLoci)
                obj.numLoci = tempNumLoci;
            else
                message = sprintf('ERROR: Number of Loci must be a positive integer\n');
            end
            %Population Size
            if (obj.matrixOn && ParameterManager.isPositiveIntegerWithin(tempPopulation, 16, obj.maxPlottingPopSize) && floor(sqrt(tempPopulation))^2 == tempPopulation) || ...
                ~obj.matrixOn && ParameterManager.isPositiveIntegerWithin(tempPopulation, 16, obj.maxPopSize);
                obj.initializeInitialCounts(obj.numLoci, tempPopulation);   
                obj.popSize = tempPopulation;
            else
            	message = sprintf('ERROR: If plotting is enabled, then population size must be a perfect square and less than %d. If plotting is not enabled, then population size must be less than %d. Population size must be at least 16.\n', obj.maxPlottingPopSize, obj.maxPopSize);

            end
            %TODO: Verify if it is better to keep this using currentType or
            %to loop over all types in this line
            %Param 1
            if ParameterManager.isNumber(tempParam1) && obj.mutating && obj.numLoci > 1
                obj.s = tempParam1;
            elseif ParameterManager.isNumberWithin(tempParam1, ...
                    obj.classConstants(obj.currentModel).ParamBounds1(1),...
                    obj.classConstants(obj.currentModel).ParamBounds1(2))
                obj.modelParameters(obj.currentModel).Param1(obj.currentType) = tempParam1;
            else
                if obj.mutating && obj.numLoci > 1
                    message = sprintf('ERROR: s must be a number\n');
                else
                    message = sprintf('ERROR: %s must be a number between %d and %d.\n', ...
                        obj.classConstants(obj.currentModel).ParamName1,...
                        obj.classConstants(obj.currentModel).ParamBounds1(1),...
                        obj.classConstants(obj.currentModel).ParamBounds1(2));
                end
            end
            
            
            %Param 2
            if ParameterManager.isNumber(tempParam2) && obj.mutating && obj.numLoci > 1
                obj.e = tempParam2;
            elseif ParameterManager.isNumberWithin(tempParam2, ...
                    obj.classConstants(obj.currentModel).ParamBounds2(1),...
                    obj.classConstants(obj.currentModel).ParamBounds2(2))
                obj.modelParameters(obj.currentModel).Param2(obj.currentType) = tempParam2;
            else
                if obj.mutating && obj.numLoci > 1
                    message = sprintf('ERROR: e must be a number\n');
                else
                    message = sprintf('ERROR: %s must be a number between %d and %d.\n', ...
                        obj.classConstants(obj.currentModel).ParamName2,...
                        obj.classConstants(obj.currentModel).ParamBounds2(1),...
                        obj.classConstants(obj.currentModel).ParamBounds2(2));
                end
            end
            
            %Ninit
            if ParameterManager.isPositiveInteger(ninitTemp)
                obj.modelParameters(obj.currentModel).Ninit(obj.currentType) = ninitTemp;
            else
                message = sprintf('ERROR: Initial Population must be a positive integer or all types.\n');
            end
            %Num Types
            if ParameterManager.isPositiveIntegerWithin(numTypesTemp, 0, obj.maxTypes)
                obj.updateNumTypes(numTypesTemp);
                obj.numTypes = numTypesTemp;
            else
                message = sprintf('ERROR: Number of types must be a positive integer.\n');
            end
            %Write any modified values back to the boxes
            obj.writeBoxes();
            
            
        end
            
        
        
        
        function updateNumTypes(obj, num)
            %Called when number of types in increased
            %Updates the parameters for all of the models
            if num ~= obj.numTypes
                if num < obj.numTypes
                    obj.handles.types_popup.String(num+1:end) = [];
                    for model = 1:length(obj.classConstants)
                        obj.modelParameters(model).Ninit(num+1:end) = [];
                        obj.modelParameters(model).Param1(num+1:end) = [];
                        obj.modelParameters(model).Param2(num+1:end) = [];
                    end
                elseif num > obj.numTypes
                    for i = obj.numTypes+1:num
                        obj.handles.types_popup.String{i} = i;
                        for model = 1:length(obj.classConstants)
                            obj.modelParameters(model).Ninit(i) = obj.modelParameters(model).Ninit_default;
                            obj.modelParameters(model).Param1(i) = obj.modelParameters(model).Param1_default;
                            obj.modelParameters(model).Param2(i) = obj.modelParameters(model).Param2_default;
                        end
                    end
                end
                %adjust birth and death rates accordingly
                if num ~= obj.numTypes
                    for i = 1:length(obj.classConstants)
                        if obj.classConstants(i).atCapacity
                            obj.modelParameters(i).Ninit = zeros(1,num);
                            for j = 1:(num-1)
                                obj.modelParameters(i).Ninit(j) = floor(obj.popSize/num);
                                obj.modelParameters(i).Ninit(j) = floor(obj.popSize/num);
                            end
                            obj.modelParameters(i).Ninit(num) = obj.popSize - sum(obj.modelParameters(i).Ninit);
                        end
                    end
                end
                if obj.numLoci == 1
                    numAlleles = num;
                else
                    numAlleles = 2;
                end
                obj.initializeMutationMatrix(numAlleles); 
            end
        end
        

        function updateNumLoci(obj, num) 
        	%Performs te necessary updates in response to the numLoci box being changed
            %Also sets the default frequencies vector
            if num ~= obj.numLoci
                if num == 1
                    numAlleles = obj.numTypes;
                else
                    numAlleles = 2;
                end
                obj.initializeMutationMatrix(numAlleles);
                obj.initializeInitialCounts(num, obj.popSize);
            end
        end

        

        function updateDefaultNinit(obj)
        	%Updates the size of the matrix/population based on the input to the
            %population box string
            for i = 1:length(obj.classConstants)
                if obj.classConstants(i).atCapacity
                    obj.modelParameters(i).Ninit = zeros(1,obj.numTypes);
                    for j = 1:(obj.numTypes-1)
                        obj.modelParameters(i).Ninit(j) = floor(obj.popSize/obj.numTypes);
                        obj.modelParameters(i).Ninit(j) = floor(obj.popSize/obj.numTypes);
                    end
                    obj.modelParameters(i).Ninit(obj.numTypes) = obj.popSize - sum(obj.modelParameters(i).Ninit);
                end
            end
        end
       
        
        function updatePopulation(obj, pop)
        	%Updates some parameters based on the input to the max
        	%population box
            obj.initializeInitialCounts(obj.numLoci, pop)
        end
       
        
        function initializeInitialCounts(obj, num, pop)
            %Initializes the initial counts vector to have all of the
            %initial organisms be of the first type
            if obj.classConstants(obj.currentModel).atCapacity
                obj.initialCounts = [pop zeros(1,2.^num - 1)];
            else
                obj.initialCounts = [1 zeros(1,2.^num - 1)];
            end
        end
        
        
        
        function initializeMutationMatrix(obj, numAlleles)
            %reinitialize mutation matrix
            obj.mutationMatrix = zeros(numAlleles);
            for i = 1:numAlleles
                for j = 1:numAlleles
                    if i == j
                        obj.mutationMatrix(i,j) = 1 - 0.01*(numAlleles-1);
                    else
                        obj.mutationMatrix(i,j) = 0.01;
                    end
                end
            end
        end
        
        
        function out = getField(obj, param)
        	%Provides an interface for accessing parameters that does not
            %require the user to know whether the numLoci > 1, or what the
            %current model is
            model = obj.currentModel;
            if strcmp(param,'numTypes')
                if ~obj.mutating || obj.numLoci == 1
                    out = obj.numTypes;
                    return;
                else
                    out = 2^obj.numLoci;
                    return;
                end
            end
            if obj.numLoci > 1 && obj.mutating
                if strcmp(param, 'Ninit')
                    assert(length(obj.initialCounts) == 2^obj.numLoci);
                    assert(sum(obj.initialCounts) <= obj.popSize);
                    out = obj.initialCounts;
%                     if ~
%                         out = [obj.modelParameters(model).Ninit_default zeros(1, 2^obj.numLoci - 1)];
%                     else
%                         tail = floor(obj.initialCounts(2:end));
%                         assert(sum(tail) <= obj.popSize);
%                         
%                         out = [obj.popSize - sum(tail) tail]; %some rounding to ensure that sum of types adds to popSize
%                     end
                elseif strcmp(param, 'Param2')
                	out = repmat(obj.modelParameters(model).Param2_default,1,2^obj.numLoci);
                elseif strcmp(param, 'Param1')
                    out = zeros(1,2^obj.numLoci);
                    for i = 1:2^obj.numLoci
                        if obj.classConstants(model).OverlappingGenerations
                            out(i) = obj.lociParam1OverlappingGenerations(i); 
                        else
                            out(i) = obj.lociParam1NonOverlappingGenerations(i);
                        end
                    end
                else
                    error('Incorrect input to getField')
                end
            else
                out = getfield(obj.modelParameters(model), param);
            end
        end

        
                
        

        function out = verifySizeOk(obj)
            %For the numLoci<1 case, verifies that the sum of the Ninits for 
            % all species is the original populaiton size 
            out = 1;
            if (obj.numLoci == 1 || ~obj.mutating)
                if obj.classConstants(obj.currentModel).atCapacity
                    if sum(obj.modelParameters(obj.currentModel).Ninit) ~= obj.popSize
                        out = 0;
                    end
                else
                    if sum(obj.modelParameters(obj.currentModel).Ninit) > obj.popSize
                        out = 0;
                    end
                end
            end
        end

        function out = getNumTypes(obj)
        	%returns the current number of types, regardless of whether
            %numLoci > 1 or not
            if obj.mutating && obj.numLoci > 1
                out = 2^obj.numLoci;
            else
                out = obj.numTypes;
            end
        end

        
        function setMatrixOn(obj, input)
        	%Sets the matrixOn variable. I might add some interesting logic
            %here
            obj.matrixOn = input;
        end
        
        function setMutationMatrix(obj, m)       
            %Sets the mutationMatrix variable.
            obj.mutationMatrix = m;
        end
        
        function setInitialCounts(obj, df)
            %sets the initialCounts variable.
            obj.initialCounts = df;
        end

        

        function out = lociParam1OverlappingGenerations(obj,num) 
        	%A function that computes the first model parameter based on the
            %number of 1s for OverlappingGenerations models (Logistic, Moren,
            %Exp, etc)
            out = 1 + obj.s*(ParameterManager.numOnes(num - 1)^(1-obj.e));
        end
        

        function out = lociParam1NonOverlappingGenerations(obj,num)
        	%A function that computes the first model parameter based on the
            %number of 1s for NonOverlappingGenerations models (Wright-Fisher)
            out = exp(obj.s*(ParameterManager.numOnes(num - 1)^(1-obj.e)));
        end
    end
        
    methods (Static) %Mostly Input Checking Functions
        function out = numOnes(x) 
            %number of one bits in input number x
            out = sum(bitget(x,1:ceil(log2(x))+1));
        end

        function out = isNumber(n)
            %verifies that input is a single real number
            out = 1;
            if ~isnumeric(n)
                out = 0;
            elseif isnan(n)
                out = 0;
            elseif ~isreal(n)
                out = 0;
            end
        end

        function out = isNumberWithin(n, lower, upper)
            %verifies that input is a single number in the range
            out = 1;
            if ~ParameterManager.isNumber(n)
                out = 0;
            elseif any(n < lower) || any(n > upper)
                out = 0;
            end 
        end

        function out = isPositiveNumber(n)
            %verifies that input is a single positive number
            out = 1;
            if ~ParameterManager.isNumber(n)
                out = 0;
            elseif any(n < 0) %negative
                out = 0;
            end 
        end

        function out = isPositiveInteger(n)
            %verifies that input is a single positive integer
            out = 1;
            if ~ParameterManager.isPositiveNumber(n)
                out = 0;
            elseif any(round(n) ~= n) %non-integer
                out = 0;
            end 
        end

       function out = isPositiveIntegerWithin(n, lower, upper)
            %verifies that input is a single number in the range
            out = 1;
            if ~ParameterManager.isPositiveInteger(n)
                out = 0;
            elseif ~ParameterManager.isNumberWithin(n, lower, upper)
                out = 0;
            end 
       end
     end
         
end
