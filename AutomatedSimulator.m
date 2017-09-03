function out = AutomatedSimulator(GridManagerClass, Ninit, p1, p2, varargin)
    % This file is an automated front end for the simulator program
    
%     Required Inputs
%     GridManagerClass - The model used for the simulation (e.g. GridManagerLogistic, GridManagerMoran) 
%     Ninit - The vector of initial populations
%     p1 - The vector of initial parameter 1s
%     p2 - The vector of initial parameter 2s
% 
%     
%     Optional Inputs
%     mutating - (whether or not mutation is enabled, Default 0)
%     numLoci - (the number of loci, Default 1)
%     mutationMatrix - (the mutation matrix, Default [0.99 0.01; 0.01 0.99])
%     recombinationNumber - (whether or not recombination is enabled,
%     Default 0)
%     matrixOn - (whether or not the petri_dish_on method is used, Default 0) 
%     totalPopSize - (the maximum population size, Default 2500)
%     spatialOn - (whether or not spatial structure is enabled, Default 1)
%     edgesOn - (whether or not edges are enabled. If 0, grid is actually a
%     torus, Default 1)
%     maxIterations - (maximum number of iterations before the
%     simulation is halted, Default 25)
%     return - (The type of data returned, Default totalCount, options are
%     totalCount, matrix, ageDist)
%     
    %Example: 
    %AutomatedSimulator('GridManagerMoran', [450 450], [1 1], [0.01 0.01], 'totalPopSize', 900)
    p = inputParser;

    addRequired(p,'GridManagerClass',@(s) exist(s, 'file'));
    addRequired(p,'Ninit',@ParameterManager.isPositiveInteger);
    addRequired(p,'p1',@ParameterManager.isNumber);
    addRequired(p,'p2',@ParameterManager.isNumber);

    addParameter(p,'mutating', 0, @(x) x == 0 || x == 1);
    addParameter(p,'mutationMatrix', [], @(M) size(M,1) == size(M,2) && ParameterManager.isPositiveNumber(M));
    addParameter(p,'numLoci', 1, @ParameterManager.isPositiveInteger);
    addParameter(p,'recombinationNumber', 0, @ParameterManager.isPositiveNumber);
    addParameter(p,'matrixOn', 0, @(x) x == 0 || x == 1);
    addParameter(p,'totalPopSize', 2500, @ParameterManager.isPositiveInteger);
    addParameter(p,'returnType', 'totalCount', @(x) any(validatestring(x,{'totalCount', 'matrix', 'ageDist'})));
    addParameter(p,'spatialOn', 0, @(x) x == 0 || x == 1);
    addParameter(p,'edgesOn', 0, @(x) x == 0 || x == 1);
    addParameter(p,'maxIterations', 25,  @ParameterManager.isPositiveInteger);

    parse(p,GridManagerClass, Ninit, p1, p2, varargin{:})
    assert(length(p.Results.Ninit) == length(p.Results.p1) && length(p.Results.p1) == length(p.Results.p2), 'ASSERTION ERROR: Lengths of Ninit, p1 and p2 must be the same!');
    
    recombining = (p.Results.recombinationNumber > 0);
    if isempty(p.Results.mutationMatrix)
        if ~p.Results.mutating || p.Results.numLoci == 1
            mutationMatrix = generateMutationMatrix(numel(p.Results.Ninit));
        else
        	mutationMatrix = generateMutationMatrix(2);
        end
    else
        mutationMatrix = p.Results.mutationMatrix;
    end
    
    if p.Results.mutating && p.Results.numLoci > 1
        assert(2^p.Results.numLoci == length(p.Results.Ninit), 'ASSERTION ERROR: Number of Loci and Number of Types is not consistent.');
        assert(all(size(mutationMatrix) == [2 2]), 'ASSERTION ERROR: Mutation matrix must be 2x2 when numLoci > 2.');
    else
    	assert(length(p.Results.Ninit) == size(mutationMatrix,1), 'ASSERTION ERROR: The sides of mutation matrix are the incorrect length.');
    end
    assert(p.Results.matrixOn || ~strcmp(p.Results.returnType, 'ageDist'), 'ASSERTION ERROR: In order to track ages, we must turn the matrix on.');
    
    %Verify that values for parameter 1 and parameter 2 are good
    if p.Results.mutating && p.Results.numLoci > 1
        assert(ParameterManager.isNumber(p.Results.p1) & ParameterManager.isNumber(p.Results.p2), 'ASSERTION ERROR: s and e must be numeric');
    else
        b1 = GUIHelper.getConstantProperty(GridManagerClass, 'ParamBounds1');
        b2 = GUIHelper.getConstantProperty(GridManagerClass, 'ParamBounds2');
        a1 = ParameterManager.isNumberWithin(p.Results.p1, b1(1), b1(2));
        a2 = ParameterManager.isNumberWithin(p.Results.p2, b2(1), b2(2));
        assert(logical(a1), 'ASSERTION ERROR: Parameter 1 (e.g. birth rate or fitness) is not within the correct range')
        assert(logical(a2), 'ASSERTION ERROR: Parameter 2 (e.g. death rate) is not within the correct range')
    end

    
    
    
    MM = MutationManager(p.Results.mutating, mutationMatrix, p.Results.numLoci, recombining, p.Results.recombinationNumber);
    
    constructor = str2func(GridManagerClass);
    gridManager = constructor(...
        p.Results.totalPopSize,...
        p.Results.Ninit,...
        MM,...
        p.Results.matrixOn,...
        p.Results.spatialOn,...
        p.Results.edgesOn,...
        p.Results.p1,...
        p.Results.p2);
    
    matCell = {};
    for iter = 1:p.Results.maxIterations
    	[c, halt] = gridManager.getNext();
        if strcmp(p.Results.returnType, 'matrix')
            matCell{iter} = gridManager.matrix;
        end
        if halt
            break;
        end
    end
    if iter >= p.Results.maxIterations
%         fprintf('Max number of iterations, %d, reached. You can set this by setting the maxIterations parameter. \n', p.Results.maxIterations);
    end
    switch p.Results.returnType
        case 'totalCount'
            out = gridManager.totalCount;
        case 'matrix'
        	out = matCell;
        case 'ageDist'
            out = gridManager.ageStructure;
    end

    function mm = generateMutationMatrix(dim)
        mm = zeros(dim);
        for i = 1:dim
            for j = 1:dim
                if i == j
                    mm(i,j) = 1 - 0.01*(dim-1);
                else
                    mm(i,j) = 0.01;
                end
            end
        end
    end

end
