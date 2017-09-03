% Multilocus stochastic mutation/selection simulator

% In unit 4 we studied the deterministic balance between recurrent
% deleterious mutation and purifying selection in a genome with an infinite
% number of loci. This piece of code will help you to explore those results
% in a finite-size genome and in the presence of random genetic drift.

% Number of generations to be simulated
maxGenerations = 1000;

% Total population size
N = 100;
% Genome size
numLoci = 10;

% mutation rate and mutation matrix: column = source and row = destination
mu = 0.025;
mutationMatrix = [1 - mu 0; ...
                  mu 1];

% selection coefficient, Moran birth and death rates.
s = 0.1;
birth = zeros(1,2^numLoci);
for i=1:2^numLoci
    birth(i) = (1-s)^MutationCount(i,numLoci);
end
death = ones(1,2^numLoci);

% The initial population configuration is roughly at deterministic equilibrium
Ninit = zeros(1,2^numLoci);
tempSum = 0;
for i = 1:2^numLoci
    Ninit(i) = ceil(N*PoissonDist(numLoci*mu/s, ...
        MutationCount(i,numLoci))/nchoosek(numLoci,MutationCount(i,numLoci)));
    tempSum = tempSum + Ninit(i);
    % if on this iteration we "overdrew" the total population, fix it.
    if tempSum > N
        Ninit(i) = Ninit(i) - (tempSum - N);
        break;
    end
end

% Call the simulator. Arguments in order:
% - model name
% - initial population configuration: wild type and then mutant
% - birth rates: wild type and then mutant
% - death rates: wild type and then mutant
% - mutations enabled flag
% - mutation matrix 
% - simulation duration
% - total population size
% - genome size
x=AutomatedSimulator('GridManagerMoran',Ninit, birth, death, ...
    'mutating',1,'mutationMatrix',mutationMatrix,'maxIterations', ...
    maxGenerations,'totalPopSize',N,'numLoci',numLoci);

% The simulator keeps track of all 2^numLoci different genotypes but since
% our fitness function is rotationally symmetric, we can pool all genotypes
% with respect to their number of deleterious mutations.
xMutationCount = zeros(numLoci+1,maxGenerations+1);
for i=1:2^numLoci
    xMutationCount(MutationCount(i,numLoci)+1,:) = ...
        xMutationCount(MutationCount(i,numLoci)+1,:) + x(i,:);
end

% ...and we also want to watch mean fitness
meanFitness = zeros(maxGenerations,1);
for i=1:maxGenerations+1
    meanFitness(i) = sum(x(:,i).*birth')/N;
end
   
% first output: kinetics of frequency for each class.
figure(1);
plot(xMutationCount'./N);
title('Kinetics of Mutant Classes');
xlabel('Generation');
ylabel('Frequency');
for i = numLoci+1:-1:1
    if sum(xMutationCount(i,:)) > 1
        i = i - 1;
        break
    end
end

% It *must* be possible to compute the right legend string using fprintf()
% but I can't quite figure out how to pass the single quotes through...
switch i
    case 0
        legend('k = 0');
    case 1
        legend('k = 0','k = 1');
    case 2
        legend('k = 0','k = 1','k = 2');
    case 3
        legend('k = 0','k = 1','k = 2','k = 3');
    case 4
        legend('k = 0','k = 1','k = 2','k = 3','k = 4');
    case 5
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5');
    case 6
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5', ...
         'k = 6');
    case 7
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5', ...
         'k = 6','k = 7');
    case 8
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5', ...
        'k = 6','k = 7','k = 8');
    case 9
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5', ...
        'k = 6','k = 7','k = 8','k = 9');
    case 10
        legend('k = 0','k = 1','k = 2','k = 3','k = 4','k = 5', ...
        'k = 6','k = 7','k = 8','k = 9','k = 10');
end

% Compute the deterministic expectation for each class. Needs to be done in
% this open-ended way (that is, we didn't pre-allocate the array) because
% we can't trivially know how long the array will be in the face of
% arbitrary mu and s values. On the other hand, we'll only collect
% counts until they drop below 1 for this N.
oldEC = 0;
i = 1;
clear ExpectedCounts
ExpectedCounts(i) = N*PoissonDist(numLoci*mu/s,0);
while oldEC < ExpectedCounts(i)
    oldEC = ExpectedCounts(i);
    i = i + 1;
    ExpectedCounts(i) = N*PoissonDist(numLoci*mu/s,i-1);
end
while ExpectedCounts(i) > 1
    i = i + 1;
    ExpectedCounts(i) = N*PoissonDist(numLoci*mu/s,i-1);
end

% second output: frequency distribution averaged over the entire run,
% together with deterministic expectations.
hold off
figure(2);
plot([0:numLoci],mean(xMutationCount,2)./N,'ko')
hold on
plot([0:i-1],ExpectedCounts./N,'rx')
hold off
title('Mean Frequency Distribution across Run');
xlabel('Mutation class (k)');
ylabel('Frequency');
legend('Observed','Deterministic expectation');

% third output: population mean fitness also averaged over the entire run.
fprintf('mean fitness = %f\n',mean(meanFitness));

% Bonus fourth output: If we're in the ratchet regime, let's see when it
% clicked
j = 1;
data = zeros(maxGenerations,1);
for i = 1:maxGenerations
    if xMutationCount(j,i) == 0
        j = j + 1;
    end
    data(i) = j-1;
end   
if data(maxGenerations) > 0
    figure(3)
    plot(data);
    title('Clicks of the Ratchet');
    xlabel('Generation');
    ylabel('Number of Mutations in the Fittest Surviving Class');
end
