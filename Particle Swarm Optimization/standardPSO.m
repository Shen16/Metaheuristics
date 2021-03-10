function [evalValues, x_tracker,stepCount, fval_tracker, pbest_tracker, gbest_tracker, back_acc_counter] = standardPSO(fhd, DIM, evalsPerDIM, bias)  %DIM= dimension, evalsPerDIM= StopStepValue or Max iterations?, fhd is rastrigin function, functionNumber= number of variable in the function==dim??
    

    evalValues = zeros(1,11);
    nextEval = 1;
    evalPoints = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    numEvals = 0;
    

    % PSO parameters -- D. Bratton and J. Kennedy, "Defining a standard for particle swarm optimization," IEEE SIS, 2007, pp. 120–127.
    popsize = 50; 
    w = 0.72984;
    c1 = 2.05 * w;  % c1 and c2 are called acceleration coefficinets or trust parameters. % c1 expresses how much confidence a particle has in itself
    c2 = 2.05 * w;  % c2 expresses how much confidence a par- ticle has in its neighbors. % w is inertia weight

    
    % Search space parameters
%    range = 5; %BBOB2009
    range = 500.25; %CEC2013  %500.25 or 5.25
    xmin = -range * ones(1,DIM); 
    xmax = range * ones(1,DIM); 
    vmin = -range * ones(1,DIM); 
    vmax = range * ones(1,DIM);

    % Random initial positions
    x = 2 * range * rand(popsize,DIM) - range;  % Modified PSO had: x = [ ones(popsize,1)*range, 2 * range * rand(popsize,DIM-1) - range ]; 1st co ordinate is fixed to range value. 
    % zero for initial velocity -- A. Engelbrecht, "Particle Swarm Optimization: Velocity Initialization," IEEE CEC, 2012, pp. 70-77.
    v = zeros(popsize,DIM);

    
    %# To Track x 
    stepCount=0; % Counter that counts the number of sets of x. So xGroupCount= 100, means 2 sets of stepCount
    xGroupCount=1;  % Counter for storing x values in rows 1:50, 51:100,...
    x_tracker= ones(20000,DIM)*50;    % initilizing 
    
    stepCount=stepCount+1;

    x_tracker(xGroupCount: (xGroupCount+ popsize-1),:)= x;    
    xGroupCount= xGroupCount+ popsize;
    
    
  
    
    % initialize personal best positions as initial positions
    pbest = x;
        % pbestCosts = feval(fhd, x);   % Uses fnd function handle with parameters x' and function Number
    pbestCosts= zeros(1, popsize);
    for i= 1:popsize
        pbestCosts(1,i) = fhd( x(i, :));   % Uses fnd function handle with parameters x' and function Number
    
    end    
    
    numEvals = numEvals + popsize;
    
    %# To count initial pbestCost and newCosts
    fval_tracker= ones(20000, popsize)* (-1);    % initilizing 
    fval_tracker(stepCount, :)= pbestCosts;      %pbestCost is a 1x50 row vector
    
    %# To count initial pbestCost only
    pbest_tracker= ones(20000, popsize)* (-1);    % initilizing 
    pbest_tracker(stepCount, :)= pbestCosts;      %pbestCost is a 1x50 row vector
    
    
    gbestCost=0; %NEW % Just Initializing to access the variable outside the loop
    
    
    %# To count initial pbestCost only
    gbest_tracker= ones(20000, 1)* (-1);   % initilizing 
    
    
    %#new for accepted (moved x) relative fitness
    %relative_fitness_acc=ones(20000, 1)* (-100); % NOT NEEEDED anymore
    
    %#new Counts the moved(accepetd) particles that are deceptive or backward steps (as in RW). (Each iteration can hace more then one particles moving since each iteration has 50 particles. We are summing all of them)
    back_acc_counter=0;
     
    % update lbest
    [lbest] = update_lbest(pbestCosts, x, popsize);

    maxevals = DIM * evalsPerDIM; 
    maxgenerations = floor(maxevals/popsize);
    for generation = 2 : maxgenerations
        
        % Update velocity 
        v = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).*(lbest-x);

        % Clamp veloctiy 
            %%% Setiing the upper and lower bounds. Anyting less than -100 is set to -100 and same for max
        oneForViolation = v < repmat(vmin,popsize,1);  %if velocity is less than -100 %oneFprViolation is a logical array, where 1 means velocity is less than -100 for a coordinate point
        v = (1-oneForViolation).*v + oneForViolation.*repmat(vmin,popsize,1); %velocities that are greater than -100 are multiplied with v (are kept) and others are just -100
        oneForViolation = v > repmat(vmax,popsize,1); 
        v = (1-oneForViolation).*v + oneForViolation.*repmat(vmax,popsize,1);   %Velocities <=100 are kep and anything greater is just 100

        
        % Update position 
        
        %#new
        init_x=x;
        
        
        x = x + v;    % Velocity added to previoys x to get new position
        
        
        % Reflect-Z for particles out of bounds -- S. Helwig, J. Branke, and S. Mostaghim, "Experimental Analysis of Bound Handling Techniques in Particle Swarm Optimization," IEEE TEC: 17(2), 2013, pp. 259-271

        % reflect lower bound
                %%% 0 or greater is the lower bound for reflection
        relectionAmount = repmat(xmin,popsize,1) - x;    % lower bound reflection= -100-x, where x is the new x
        oneForNeedReflection = relectionAmount > zeros(popsize,DIM);
        relectionAmount = (1-oneForNeedReflection).*zeros(popsize,DIM) + oneForNeedReflection.*relectionAmount;  % The reflection amount that are less than or equal to 0 are made 0 and others are kept 
        % clampfirst
        x = (1-oneForNeedReflection).*x + oneForNeedReflection.*repmat(xmin,popsize,1); % keep the x for the reflection that are less than 0 and make others -100
        
        % then reflect
        x = x + relectionAmount;  % x+0 for oneForNeedreflection <0 and relectionamount+ (-100) for reflection amount >0
        

        % set velocity for reflected particles to zero
        v = (1-oneForNeedReflection).*v + oneForNeedReflection.*zeros(popsize,DIM);  %Particles with reflectionamount> 0 are called reflected particles and their velocities is set to 0
        
        % reflect upper bound
        relectionAmount = repmat(xmax,popsize,1) - x;   % upper bound reflection= 100-x, where x is the new x after lower bound
        oneForNeedReflection = relectionAmount < zeros(popsize,DIM);
        relectionAmount = (1-oneForNeedReflection).*zeros(popsize,DIM) + oneForNeedReflection.*relectionAmount; % reflectionamount that are greater than equal to zero are made 0 and others are kept
        % clampfirst
        x = (1-oneForNeedReflection).*x + oneForNeedReflection.*repmat(xmax,popsize,1); % keep the x for the reflection that are greater than 0 and make others 100
        % then reflect
        x = x + relectionAmount; % x+0 for oneForNeedreflection >0 and relectionamount+ (100) for reflection amount <0

        % set velocity for reflected particles to zero
        v = (1-oneForNeedReflection).*v + oneForNeedReflection.*zeros(popsize,DIM); %Particles with reflectionamount< 0 are called reflected particles and their velocities is set to 0
        
        
        
        new_x= x;
        basinFitness= @(x) sum((bias*round(abs(x)))) ;
        
        for i= 1:popsize
            if (fhd(new_x(i, :))< fhd(init_x(i, :)) && basinFitness(new_x(i,:))> basinFitness(init_x(i,:)))
                back_acc_counter= back_acc_counter+1;  
            end
        end
        
        
        
        
        %{
        
        %#new %NOT NEEDED anymore!
        %to find if a particle moved (out of 50 particles)
        new_x= x;
        energy= ones(1,popsize)*-100;
        for i= 1:popsize
            if (new_x(i)~= init_x(i))
                energy(i)= fhd( new_x(i, :));  %50 rows or partivles's energy/fitness so 50 enegergies or less if a particle did not move
            end
        end 
        
        energy_valid= energy(energy~=-100); %getting energies for particles that moved
        min_energy= min(energy_valid); %grtting the min energy 
        
        %Finding the particle with minimum energy
        min_x=zeros(1,DIM);
        for i= 1:popsize
            if (fhd( new_x(i, :))== min_energy)
                min_x(1,:)= new_x(i, :);  %50 rows or partivles's energy so 50 enegergies or less if a particle did not move
            end
        end 
        
        %Getting rid of the fitness of the basin
        basinFitness= @(x) sum((bias*round(abs(x)))) ;
        relative_fitness_acc(stepCount)= fhd( min_x(1,:))- basinFitness(min_x(1,:));
        %end of #new code
        
        %}
        
        
        
        
        x_tracker(xGroupCount: (xGroupCount+ popsize-1),:)= x;    
        xGroupCount= xGroupCount+ popsize;
        stepCount=stepCount+1;
        
        % Update pbest, lbest, and gbest
                %newCosts = feval(fhd, x');
        newCosts= zeros(1, popsize);
        for i= 1:popsize
            newCosts(1,i) = fhd( x(i, :));   % Uses fnd function handle with parameters x' and function Number
        end   
        
        numEvals = numEvals + popsize;

        fval_tracker(stepCount, :)= newCosts;  
        
        
        for index = 1:popsize
            if newCosts(index) < pbestCosts(index)
                pbest(index,:) = x(index,:);
                pbestCosts(index) = newCosts(index);
            end
        end
        pbest_tracker(stepCount, :)= pbestCosts; 
        
        [lbest] = update_lbest(pbestCosts, pbest, popsize);
        
        
        
        
        %{
        %To compute relative frequency
        if (gbestCost~=min(pbestCosts))
            updated_gbest_tracker(stepCount-1)= gbestCost; %only when gbest is different from previous
        end
        
        %}
        
        
        gbestCost = min(pbestCosts);
        
        gbest_tracker(stepCount-1)= gbestCost;

        if numEvals >= evalPoints(nextEval) * maxevals
            evalValues(1, nextEval) = gbestCost;
            nextEval = nextEval + 1;
        end
    end
    evalValues(1, 11) = gbestCost;
end


% Function to update lbest
function [lbest] = update_lbest(costX, x, popsize)
    %particle 1 is neighbours with particle n=popsize
    sm(1, 1)= costX(1, popsize);
    sm(1, 2:3)= costX(1, 1:2);
    [cost, index] = min(sm);
    
    %%%This part is to get the first row of lbest values
    if index==1
        lbest(1, :) = x(popsize, :);
    else
        lbest(1, :) = x(index-1, :);
    end

    %%%This part is to get the next 2nd row to 49th row of lbest values
    for i = 2:popsize-1
        sm(1, 1:3)= costX(1, i-1:i+1);
        [cost, index] = min(sm);
        lbest(i, :) = x(i+index-2, :);
    end

    %%%This part is to get the 50th row of lbest values
    % particle n=popsize is neighbours with particle 1
    sm(1, 1:2)= costX(1, popsize-1:popsize);
    sm(1, 3)= costX(1, 1);
    [cost, index] = min(sm);
    if index==3
        lbest(popsize, :) = x(1, :);
    else
        lbest(popsize, :) = x(popsize-2+index, :);
    end    
end
