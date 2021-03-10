%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   anneal.m  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Usage:
%     [MINIMUM,FVAL] = ANNEAL(LOSS,NEWSOL,[OPTIONS]);
%          MINIMUM is the solution which generated the smallest encountered
%          value when input into LOSS.
%          FVAL is the value of the LOSS function evaluated at MINIMUM.



%diary 1D_SA_Rastr
%%%%%%% 
% Rastrigin function handle with anneal


tic; %start time

d=10;  %DIM
n=30;   % Number of simulations


basin_height=600;
bias=3;
evalsPerDIM=10000;
popSize=50;

fval= zeros(n,1);  % To record the final evalValye after each simulation

GroupCount= 0;

fbest_index= 100; %Just intiaizing
gbest_x=  zeros(n,d);


%#new
mean_rf_acc= zeros(n,1); %mean values of rf_acc for each n cycles

back_acc_count= zeros(n,1); %Storing counts of accepted backward (or deceptive) steps for n iterations


for i=1:n
    GroupCount= GroupCount +popSize;
    
    rastrigin = @(x) sum( (bias*round(abs(x))) - (basin_height/2)/d* (cos(2*pi*x)-1)); % x is a row vector of length d  %% make cos part 0 to 30 again (decides depth)  % global (0,0)
    
    %Intial position parent
    %{
    a = -5.12;
    b = 5.12;
    parent = (b-a).*rand(1,d) + a; % Generates a random intial position of d dimension, within +5.12 and -5.12 range
    %}
    
    %parent= ones(1,d)* 500.25; % Starts from 500.25 (with InitS=5) or 5.25 (with InitS=1)
    
    
    
    [evalValues, xT, stepCount, fT, pT, gT, back_acc] = standardPSO(rastrigin, d, evalsPerDIM, bias); % evals per dim 3000 before
    
    
    %# To generate Distance from origin histogram
    
    arrayLength_f = length(fT(fT(:,1)~= -1));
    f= fT(1:arrayLength_f, :);
    fbest_val= min( f(end, :) );
    for l=1:popSize
        if f(end,l)==fbest_val
            fbest_index= l;
        end
    end
    
    
    arrayLength = length(xT(xT(:,1)~= 50));
    x_all= xT(1:arrayLength, :);
    
    no_of_sets= arrayLength/popSize;
    for  z= (arrayLength - (popSize-1)) : arrayLength
        row= arrayLength- (popSize-1) + (fbest_index-1);  %index is the f_valtracker index that is minimum
        gbest_x(i,:) =x_all(row, : ) ;
    end
    
    %{
    dist_individual= zeros(popsize, 1);
    origin=zeros(1,d);
    for j=1:popsize
        dist_individual(j)= sqrt(sum((x_all(end-(popsize-j))-origin).^2));
    end
    
    dist_frm_origin (i)= mean(dist_individual);
    %}
    
    %#To generate the f(x) histogram
    fval(i)=evalValues(end,end);  %Last generation of each evalValues (which has 11 generations) is the final fval
   
    
    %{
    %NOT NEEDED anymore
    %#new
    %Relative Fitness for accepted particles (choosing the particle with minimum energy out of 50 particles (==pbest if all particles move??)). Accepted means that the particle was moved (or velocity~=0) 
    rf_acc_valid=rf_acc(rf_acc~=-100);
    mean_rf_acc(i)=mean(rf_acc_valid);
    %}
    
    %#new
    %Counts of accepted backward (or deceptive) steps
    back_acc_count(i)= back_acc;
    
end


time= toc; %end timer



name= sprintf('%dD_PSO_Rastr', d);
filename= name + ".mat";
save(filename, 'd', 'basin_height', 'bias', 'evalsPerDIM', 'popSize', 'n', 'fval', 'GroupCount','fbest_index', 'gbest_x', 'rastrigin', 'evalValues', 'xT', 'stepCount', 'fT', 'pT', 'gT', 'time', 'back_acc_count' )




%{
%#
%Average Relative frequency of accepted particles 
mean(mean_rf_acc)

Dimension = [1;2;3;5;10;20;30];
Mean = [2.0925;8.1591;13.0579;19.2221;20.8217;21.4194;21.8569];


PSO = table(Dimension, Mean)
%}

load 

%Mean fitness table of PSO 

Dimension = [1;2;3;5;10;20;30];
Mean_fitness = [1.9251e-13; 0.22024; 3.689; 14.674; 41.327; 91.768; 133.72];


PSO_mean_fitness = table(Dimension, Mean_fitness)


plot(Dimension, Mean_fitness)


title("Mean Fitness")
xlabel('Dimension','fontSize',12);
ylabel('Fitness','fontSize',12);











