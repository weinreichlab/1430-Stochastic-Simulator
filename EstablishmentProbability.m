% This program uses the AutomatedSimulator() interface to our stochastic
% simulator to perform replicate realizations of an exponentially growing
% population. One purpose might be to explore the relationship between
% analytic predictions for establishment probability developed in class and
% simulation results. But you are strongly encouraged to let your
% imagination take you elsewhere...! 
%
% September 13, 2015

% This determines the number of replicate simulations performed
Replicates = 1000;
% This is the max number of generations the simulator will run
maxIterations = 200;
% This is the maximum possible value of population size
Nmax = 1000;

% This is the number of individuals with which the simulation starts
InitialCount = 1;
% These are the instantaneous per-capita probabilities of reproduction and
% death
b = 1;
d = 0.75;

% This matrix will receive our results: the number of organisms at each time. It has one row for each replicate,
% and one column for each generation. zeros() intializes a matrix of given
% dimensions with zero in each cell. 
MyData = zeros(Replicates,maxIterations);

% This counter will be used to keep track of the number of
% replicates in which the population got established.
Established = 0;

% Here's the main loop: our index i will count up from 1 to Replicates
for i=1:Replicates
    % Here's where we do the work: call the simulator and capture the
    % results in a row vector. (Matlab doesn't insist that
    % variables be created or initialized before the first time they're
    % used.) The length of the row vector will vary from one
    % replicate simulation to another.  If you want to o explore the syntax and
    % capabilities of AutomatedSimulator(), open AutomatedSimulator.m and
    % take a look at the comments near the top of that file.
    ThisRun = AutomatedSimulator('GridManagerExp',[InitialCount],[b],[d],'totalPopSize',Nmax,'maxIterations',maxIterations);

    % Now we want to copy ThisRun into MyData.  
    % The realization will have turned out in one of three ways: either
    % the lineage went extinct, the lineage reached maxPopSize, or the
    % lineage remained at some intermediate size (unlikely).  The folowing
    % if/else statements handle those different eventualities.
    
    if ThisRun(end) == 0     % population went extinct
        MyData(i, 1:length(ThisRun)) = ThisRun; % remaining columns remain zero
    elseif ThisRun(end) == Nmax   % population got up to maximum allowed value
        MyData(i,1:length(ThisRun)) = ThisRun;
        MyData(i, 1+ length(ThisRun): end) = ThisRun(end);
        Established = Established + 1;
    else         % population remained at intermediate size (unlikely)
        MyData(i,1:length(ThisRun)) = ThisRun;
        % students: what is the logic behind the following if statement?
        if (d/b)^MyData(end) < 0.5
            Established = Established + 1;
        end 
    end
end

% Let's plot #organisms vs. t for those trials that didn't go extinct. 
t = 0:maxIterations - 1;
notExtinct= find( MyData( :, end) > 0);  % trials where final # orgs > 0
Extinct = find( MyData( :, end) == 0);
if size(notExtinct,1) > 0
    figure (1)
    semilogy( t, MyData(notExtinct, : ) )
    xlabel('time (generations)')
    ylabel('# organisms')
    % you might want to change the x-axis by uncommenting the following line
    % xlim([0,30])
    figure (2)
    plot( t, MyData(Extinct, : ) )
    xlabel('time (generations)')
    ylabel('# organisms')
    x = [0:max(MyData(:,10))];
    for i=1:10
        figure (i+2)
        histogram(MyData(:,i),x);
         ylim([0,Replicates])
    end

    sprintf('%d out of %d replicate realizations established.',Established,Replicates)
else
    sprintf('no replicates established.')
end
