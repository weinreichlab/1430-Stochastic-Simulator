classdef (Abstract) GridManagerAbstract < handle
%This class is the parent class of the 4 grid managers. It contains methods
%that are used by all of its child classes, as well as the instance
%variables used across all classes. 

%The GridManager class stores the matrix/counts of each species, and
%contains a method getNext that advances the matrix/counts to the next
%generation, based on the model that the GridManager's implementation is based
%on
    properties (Abstract, Constant)  %The tag properties, these characterize the class itself
        Name; %The name of the model
        OverlappingGenerations; %Whether the generations overlap. This is 1 for Moran, 0 for Wright-Fisher
        ParamName1; %Name of first parameter, e.g. birth_rate
        ParamName2; %Name of second parameter, e.g. birth_rate
        atCapacity; %Whether or not the model requires that the total population is always constant and at capacity. 1 for moran, 0 for logistic
        plottingEnabled; %Whether or not the model supports a matrixOn (show petri dish) mode
    end 
    
    properties (SetAccess = private)
        maxSize; %The maximum possible population size
        mutationManager; %A pointer to the mutation manager class that handles mutating each type to other types
        matrixOn; %Determines whether the matrix model or vector model is used
        spatialOn; %Determines whether the spatial structure is considered. Functionality is determined by child classes.
        edgesOn; %Determines whether the edges are considered when searching for the closest square of a type. If 0, the matrix is actually a torus
        Param1; %The first parameter for the model. Typically birth rate or fitness.
        Param2; %The second parameter for the model. Typically death rate. Not all models have 2 parameters. 
        matrix; %Stores the location of each organism
        oldMatrix; %Stores the matrix from the previous iteration
        ageMatrix; %stores how old each organism in the matrix is 
        ageStructure; %cell array, each cell is a timestep storing a vector that stores the frequency of each age
        saveData; %The data that is exported to a .mat file when the user presses save
        timestep; %The number of steps since the beginning of the simulation
        numTypes; %The number of different types in the simulation
        colors; %The vector of each type's color

    end
    
    properties (SetAccess = protected)
        %The Population Parameters. Updated in updateParams. 
        totalCount; %Stores the counts of each type in each generation. matrix of dimensions num types x timestep
        percentCount; %Stores the percent counts of each type in each generation. matrix of dimensions num types x timestep
        overallMeanFitness; %Stores the overall fitness of each type in each generation. vector of length timestep
    end
    
    methods (Access = public)
        %The constructor method. This function initializes the matrix and
        %the plotting parameters, as well as other useful variables like
        %the color variable. The final 2 inputs are the Param1 and the
        %Param2 inputs
        function obj = GridManagerAbstract(maxSize, Ninit, mutationManager, matrixOn, spatialOn, edgesOn, p1, p2)
            assert(~matrixOn || floor(sqrt(maxSize))^2 == maxSize)
            assert((obj.atCapacity && sum(Ninit)==maxSize) || (~obj.atCapacity && sum(Ninit) <= maxSize), 'ASSERTION ERROR: Incorrect initial populations');
            obj.Param1 = p1';
            obj.Param2 = p2';
            obj.spatialOn = spatialOn;
            obj.matrixOn = matrixOn;
            obj.edgesOn = edgesOn;
            if obj.matrixOn
                obj.matrix = zeros(sqrt(maxSize));
                obj.ageMatrix = zeros(sqrt(maxSize)) - 1;
            else
                obj.matrix = [];
                obj.ageMatrix = [];
            end
            obj.oldMatrix = obj.matrix;
            obj.maxSize = maxSize;
            obj.timestep = 1;
            obj.numTypes = length(Ninit);
            obj.totalCount = Ninit';
            
            
            if obj.matrixOn
                if (sum(Ninit) == obj.maxSize) || ~spatialOn %static population models, random placement
                    r = randperm(numel(obj.matrix));
                    i = 1;
                    for type = 1:length(Ninit)
                        n = Ninit(type);
                        obj.matrix(r(i:n+i-1)) =  type;
                        i = i + n;
                    end
                else %non-static models, place founding cells in center
                    origVec = [];
                    for type = 1:length(Ninit)
                        origVec = [origVec repmat(type,1,Ninit(type))];
                    end
                    origVec = origVec(randperm(length(origVec)));
                    ind = obj.getCenter();
                    for i = 1:length(origVec)
                        obj.matrix(ind) = origVec(i);
                        [a, b] = ind2sub(size(obj.matrix), ind);
                        ind = obj.getNearestFree(a,b);
                    end
                end
            end
           

            %when there are more than 10 types, colors become randomized
            obj.colors = [1 0 0; ...
                0 1 0; ...
                0 0 1; ...
                1 1 0; ...
                1 0 1;...
                0 1 1; ...
                1 0.5 0.5; ...
                0.5 1 0.5; ...
                1 0.5 0.5;...
                0.5 0.5 0.5; ...
                0.25 0.25 1; ...
                1 .25 .25;...
                .25 1 .25;
                .5 .25 0; ...
                .25 0 .5; ...
                0 .25 .5; ...
                .15 .15 .15;];
            obj.percentCount = [];
            obj.overallMeanFitness = [];
            obj.ageStructure = {};
            obj.mutationManager = mutationManager;
            obj.saveData = struct('Param1', p1, 'Param2', p2, 'matrix', obj.matrix);
            obj.updateParams();
        end
      
                
        
        %This method updates obj.totalCount for the new timestep, and, if
        %matrixOn is enabled, also updates the GridManager's petri dish
        %changed - entries in matrix that have changed
        %h - whether or not we should halt
        function [changed, h] = getNext(obj)
            assert(min(obj.totalCount(:, obj.timestep)) >= 0);
            obj.totalCount(:, obj.timestep + 1) = obj.totalCount(:, obj.timestep);
            obj.timestep = obj.timestep + 1;
            obj.getNextGeneration();
            %             obj.mutationManager.mutate(obj);
            %             obj.mutationManager.recombination(obj);
            if ~obj.matrixOn
                changed = [];
            else
                changed = find(obj.oldMatrix ~= obj.matrix);
                obj.oldMatrix = obj.matrix;
            end
            obj.updateParams();
            obj.saveData.totalCount = obj.totalCount;
            obj.saveData.totalCount = obj.totalCount;
            obj.saveData.ageStructure = obj.ageStructure;
            obj.saveData.matrix = cat(3, obj.saveData.matrix, obj.matrix);
            h = max(obj.totalCount(:, obj.timestep))>=obj.maxSize;
            h = h && ~obj.mutationManager.mutating || (sum(obj.totalCount(:, obj.timestep)) == 0);
        end
        
        
        %The default implementation for getNextGeneration for an
        %OverlappingGeneration models
        function getNextGeneration(obj)
            assert(obj.OverlappingGenerations == 1, 'ERROR: getNext not overriden by a non-Overlapping Generations class');
            for i = 1:sum(obj.totalCount(:, obj.timestep))
                obj.reproductiveEvent();
                if max(obj.totalCount(:, obj.timestep)) == 0
                    break;
                end
            end
            %then, include all computation updates
        end  
        
        %A method that should be left unimplemented in any
        %non-OverlappingGenerations class, and should be implemented in any
        %OverlappingGenerations class
        function reproductiveEvent(obj)
            assert(obj.OverlappingGenerations == 0, 'ERROR: reproductiveEvent called in an OverlappingGenerations class');
        end


        
        %Returns a free square in the matrix
        function ind = getFree(obj)
            assert(obj.matrixOn == 1)
            free = find(obj.matrix == 0);
            if isempty(free)
                ind = randi(obj.maxSize);
            else       
                ind = free(randi(length(free)));
            end
        end

        %Returns the nearest square in the matrix of type t        
        function ind = getNearestOfType(obj, x, y, t)
            assert(obj.matrixOn == 1);
            indices = find(obj.matrix == t);
            base = sub2ind(size(obj.matrix), x, y);
            indices = indices(indices ~= base);
            if isempty(indices)
                ind = 0;
                return;
            end
            dists = matrixDistance(obj, repmat(base,1,length(indices))', indices);
            [~,i] = min(dists);
            ind = indices(i);
        end
        
        %Returns the manhattan distance between 2 cells in the matrix. Uses wrapping
        %if ~obj.edges
        function d = matrixDistance(obj, base, ind)
            assert(obj.matrixOn == 1);
            [a_1, b_1] = ind2sub(size(obj.matrix), base);
            [a_2, b_2] = ind2sub(size(obj.matrix), ind);
            wd = abs(a_1 - a_2);
            hd = abs(b_1 - b_2);
            if ~obj.edgesOn %wrappping
                wd = min(wd, size(obj.matrix,1) - wd);
                hd = min(hd, size(obj.matrix,2) - hd);
            end
            d = hd + wd;
        end
        

        
        %Returns the nearest free square in the matrix
        function ind = getNearestFree(obj, i, j)
            assert(obj.matrixOn == 1);
            ind = getNearestOfType(obj, i, j, 0);
        end
        
        %Returns the non-diagonal neighboring squares of i,j
        % indices(1,:) = row
        % indices(2,:) = column
        %     4
        %   1 x 3
        %     2
        function indices = getNeighbors(obj, i, j)
            assert(obj.matrixOn == 1);
            s = sqrt(obj.maxSize);
            %unchanged
            indices(1,1) = i;
            indices(1,3) = i;
            indices(2,2) = j;
            indices(2,4) = j;
            %changed
            indices(2,1) = j - 1;
            indices(2,3) = j + 1;
            indices(1,2) = i + 1;
            indices(1,4) = i - 1;
            if ~obj.edgesOn
            	indices(2,1) = mod(indices(2,1) - 1, s) + 1;
                indices(2,3) = mod(indices(2,3) - 1, s) + 1;
                indices(1,2) = mod(indices(1,2) - 1, s) + 1;
                indices(1,4) = mod(indices(1,4) - 1, s) + 1;
            else
                remove = [];
                if indices(2,1) <= 0 || indices(2,1) > s
                    remove = [remove 1];
                end
                if indices(2,3) <= 0 || indices(2,3) > s
                    remove = [remove 3];
                end
                if indices(1,2) <= 0 || indices(1,2) > s
                    remove = [remove 2];
                end
                if indices(1,4) <= 0 || indices(1,4) > s
                    remove = [remove 4];
                end
                indices(:,remove) = [];
            end
            assert(length(indices) >= 2 && length(indices) <= 4, 'ERROR: Number of Neighbors is incorrect');
        end
        
        
        %Get an occupied neighbor of the input cell.
        function RowCol = getNeighborOccupied(obj, a, b)
            neighbors = obj.getNeighbors(a, b);
            neighbors = neighbors(:, randperm(size(neighbors, 2))); %prevent preferential treatment
            weights = zeros(1, length(neighbors));
            index = obj.weightedSelection(weights);
            RowCol = [neighbors(1,index), neighbors(2,index)];
        end
        
        
        %Get the neighbors of the input cell. If any neighbor is
        %free, select it. Otherwise, randomly select a
        %neighbor weighted by typeWeighting. 0 cells are always chosen if
        %possible
        %typeWeighting - a vector of length numTypes that assigns a weight
        %to the chance of killing each type
        function RowCol = getNeighborWeighted(obj, a, b, typeWeighting)
            neighbors = obj.getNeighbors(a, b);
            neighbors = neighbors(:, randperm(size(neighbors, 2))); %prevent preferential treatment
            weights = zeros(1, length(neighbors));
            for w = 1:length(weights)
                t = obj.matrix(neighbors(1,w), neighbors(2,w));
                if t == 0
                    RowCol = [neighbors(1,w), neighbors(2,w)];
                    return;
                else
                    weights(w) = typeWeighting(t);
                end
            end
            if all(weights == 0)
                index = randi(length(weights));
            else
                index = obj.weightedSelection(weights);
            end
            RowCol = [neighbors(1,index), neighbors(2,index)];
        end

        %Gets the center cell in the matrix
        function ind = getCenter(obj)
        	assert(obj.matrixOn == 1);
            a = ceil(size(obj.matrix, 1)/2);
            ind = sub2ind(size(obj.matrix),a,a);
        end
        
        %Gets the ith color from the color matrix
        function c = getColor(obj,i)
            c = obj.colors(i,:);
        end

        %returns 1 if there is only one species
        function h = isHomogenous(obj)
            h = (max(obj.totalCount(:, obj.timestep)) == obj.maxSize);
        end
        
        %Returns a random cell of the chosen type
        function out = getRandomOfType(obj, type)
            assert(obj.matrixOn == 1)
            ofType = find(obj.matrix == type);
            if isempty(ofType)
                out = -1;
            else
            	out = ofType(randi(numel(ofType)));
            end
        end
                
        %Returns an index in the vector vec, weighted by the contents of
        %vec
        function [ind, num] = weightedSelection(obj, vec)
            %TODO: Make sure that weights are never negative or all zero
            %do weighted selection
            num = rand()*sum(vec);
            ind = 0;
            while num > 0
                ind = ind + 1;
                num = num - vec(ind);
            end
        end
        
        %Returns a random type from among the valid types for this matrix
        function type = getRandomType(obj)
            type = randi(obj.num_types);
        end
       
        %Externally facing Birthing function
        % - Adds an organism to matrix
        % - Updates ageMatrix
        % - Performs mutation
        % - Updates totalCount for current timestep
        % - If matrixOn and the place that we are birthing to is occupied,
        %   then this function calls kill
        % ind - index in matrix to place the new cell
        % type - new type
        function birth(obj, ind, type)
            assert((ind == 0 && ~obj.matrixOn) || (ind > 0 && obj.matrixOn), sprintf('index: %d type: %d', ind, type));
            assert((type > 0) && (type <= obj.numTypes) && (round(type) == type), 'ERROR: New type must be an integer between 1 and numTypes');
            obj.totalCount(type, obj.timestep) = obj.totalCount(type,obj.timestep) + 1;
            if obj.matrixOn
                %If cell is occupied, kill whats currently in it
                if obj.matrix(ind) ~= 0
                    obj.kill(ind, obj.matrix(ind));
                end
                obj.matrix(ind) = type;
                obj.ageMatrix(ind) = 0;
            end
            newType = obj.mutationManager.atomicMutation(obj, ind, type);
            obj.mutationManager.atomicRecombination(obj, ind, newType);
        end
        
        % Externally facing Killing function
        % - Removes an organism from matrix
        % - Updates ageMatrix
        % - Updates totalCount for current timestep
        % ind - index in matrix of dead cell (0 if ~matrixOn) 
        % type - dead type
        function kill(obj, ind, type)
            assert((ind == 0 && ~obj.matrixOn) || (ind > 0 && obj.matrixOn && obj.matrix(ind) == type));
            assert((type > 0) && (type <= obj.numTypes) && (round(type) == type), 'ERROR: Kill type must be an integer between 1 and numTypes');
            assert(obj.totalCount(type, obj.timestep) > 0, 'ERROR: Total count matrix is already 0 at the kill type');
            obj.totalCount(type, obj.timestep) = obj.totalCount(type,obj.timestep) - 1;
            if obj.matrixOn
                obj.matrix(ind) = 0;
                obj.ageMatrix(ind) = -1;
            end
        end
        
        
        %Externally facing mutating function - called by MutationManager
        %functions
        % - Changes cell in matrix
        % - Updates totalCount for current timestep
        % ind - index in matrix of mutating cell (0 if ~matrixOn) 
        % oldType - current type of mutating cell
        % newType - new type of mutating cell
        function mutate(obj, ind, oldType, newType)
            assert((ind == 0 && ~obj.matrixOn) || (ind > 0 && obj.matrixOn), 'ERROR: ind not 0 on matrix off or  ind 0 on matrix on');
            assert(~obj.matrixOn || obj.matrix(ind) == oldType, 'ERROR: ind and oldType dont match up');
            assert((oldType > 0) && (oldType <= obj.numTypes) && (round(oldType) == oldType), 'ERROR: Old type must be an integer between 1 and numTypes');
            assert((newType > 0) && (newType <= obj.numTypes) && (round(newType) == newType), 'ERROR: New type must be an integer between 1 and numTypes');
            assert(obj.totalCount(oldType, obj.timestep) > 0, 'ERROR: There are no cells of the oldType');
            obj.totalCount(oldType, obj.timestep) = obj.totalCount(oldType,obj.timestep) - 1;
            obj.totalCount(newType, obj.timestep) = obj.totalCount(newType,obj.timestep) + 1;
            if obj.matrixOn
                obj.matrix(ind) = newType;
            end
        end
        
        % Resets the matrix and total count to the input
        % Resets the age matrixc
        % Mutates the matrix
        function newMatVec(obj, newMat, newVec)
            assert(sum(newVec) <= obj.maxSize, 'ERROR: The sum of newVecs elements is too large');
            assert(all(newVec >= 0), 'ERROR: newVec has negative elements');
            if obj.matrixOn
                assert(all(size(obj.matrix) == size(newMat)), 'ERROR: Dimensions of newMat are wrong');
                assert(all(newMat(:) >= 0) && all(newMat(:) <= obj.numTypes) && all(round(newMat(:)) == newMat(:)), 'ERROR: newMat has invalid elements');
                obj.matrix = newMat;
                obj.ageMatrix = zeros(sqrt(obj.maxSize));
                obj.ageMatrix(newMat == 0) = -1;
            end
            obj.totalCount(:, obj.timestep) = newVec;
            obj.mutationManager.mutateGeneration(obj);
            obj.mutationManager.recombineGeneration(obj);
        end
                    
        %Updates the totalCount, the percentCount, the overallMeanFitness and
        %the ageMatrix/ageStructure properties
        function updateParams(obj)
            meanFitness = zeros(1,obj.numTypes);
            for i = 1:obj.numTypes
                if sum(obj.totalCount(:, obj.timestep)) == 0
                    obj.percentCount(i, obj.timestep) = 0;
                else
                    obj.percentCount(i, obj.timestep) = obj.totalCount(i, obj.timestep)./sum(obj.totalCount(:, obj.timestep));
                end
            end
            obj.overallMeanFitness(obj.timestep) = dot(obj.Param1, obj.percentCount(:, obj.timestep));
            %dot(meanFitness, obj.totalCount(:,obj.timestep))
            if obj.matrixOn
                obj.ageMatrix(obj.ageMatrix ~= -1) = obj.ageMatrix(obj.ageMatrix ~= -1) + 1;
                maxAge = max(1, max(obj.ageMatrix(:)) + 1);
                obj.ageStructure{obj.timestep} = zeros(obj.numTypes, maxAge);
                for i = 1:obj.numTypes
                    ages = obj.ageMatrix(obj.ageMatrix ~= -1 & obj.matrix == i);
                    ageHist = hist(ages, max(ages))./length(ages);
                    obj.ageStructure{obj.timestep}(i,:) = [ageHist zeros(1,maxAge - length(ageHist))];
                end
            end            
        end

    end
    
    
    

    
    
end