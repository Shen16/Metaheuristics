function [minimum, fval, rejs, xTracker, fvalTracker, All_energies, acc_energies, acc_it , last_BAccR, sR, sBlk,sBlue, Max_S, rf_accepted, back_acc_count] = anneal(loss, parent, options, d , bias) 
% ANNEAL  Minimizes a function with the method of simulated annealing
% (Kirkpatrick et al., 1983)
%
%  ANNEAL takes three input parameters, in this order:
%
%  LOSS is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARENT is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing. If no
%  OPTIONS structure is provided, ANNEAL uses a default structure. OPTIONS
%  can contain any or all of the following fields (missing fields are
%  filled with default values):
%
%       Verbosity: Controls output to the screen.
%                  0 suppresses all output
%                  1 gives final report only [default]
%                  2 gives temperature changes and final report 
%       Generator: Generates a new solution from an old one.
%                  Any function handle that takes a solution as input and
%                  gives a valid solution (i.e. some point in the solution
%                  space) as output.
%                  The default function generates a row vector which
%                  slightly differs from the input vector in one element:
%                  @(x) (x+(randperm(length(x))==length(x))*randn/100)
%                  Other examples of possible solution generators:
%                  @(x) (rand(3,1)): Picks a random point in the unit cube
%                  @(x) (ceil([9 5].*rand(2,1))): Picks a point in a 9-by-5
%                                                 discrete grid
%                  Note that if you use the default generator, ANNEAL only
%                  works on row vectors. For loss functions that operate on
%                  column vectors, use this generator instead of the
%                  default:
%                  @(x) (x(:)'+(randperm(length(x))==length(x))*randn/100)'
%        InitTemp: The initial temperature, can be any positive number.
%                  Default is 1.
%        StopTemp: Temperature at which to stop, can be any positive number
%                  smaller than InitTemp. 
%                  Default is 1e-8.
%         StopVal: Value at which to stop immediately, can be any output of
%                  LOSS that is sufficiently low for you.
%                  Default is -Inf.
%       CoolSched: Generates a new temperature from the previous one.
%                  Any function handle that takes a scalar as input and
%                  returns a smaller but positive scalar as output. 
%                  Default is @(T) (.8*T)
%      MaxConsRej: Maximum number of consecutive rejections, can be any
%                  positive number.
%                  Default is 1000.
%        MaxTries: Maximum number of tries within one temperature, can be
%                  any positive number.
%                  Default is 300.
%      MaxSuccess: Maximum number of successes within one temperature, can
%                  be any positive number.
%                  Default is 20.
%
%
%  Usage:
%     [MINIMUM,FVAL] = ANNEAL(LOSS,NEWSOL,[OPTIONS]);
%          MINIMUM is the solution which generated the smallest encountered
%          value when input into LOSS.
%          FVAL is the value of the LOSS function evaluated at MINIMUM.
%     OPTIONS = ANNEAL();
%          OPTIONS is the default options structure.
%
%
%  Example:
%     The so-called "six-hump camelback" function has several local minima
%     in the range -3<=x<=3 and -2<=y<=2. It has two global minima, namely
%     f(-0.0898,0.7126) = f(0.0898,-0.7126) = -1.0316. We can define and
%     minimise it as follows:
%          camel = @(x,y)(4-2.1*x.^2+x.^4/3).*x.^2+x.*y+4*(y.^2-1).*y.^2;
%          loss = @(p)camel(p(1),p(2));
%          [x f] = ANNEAL(loss,[0 0])
%     We get output:
%               Initial temperature:     	1
%               Final temperature:       	3.21388e-007
%               Consecutive rejections:  	1027
%               Number of function calls:	6220
%               Total final loss:        	-1.03163
%               x =
%                  -0.0899    0.7127
%               f =
%                  -1.0316
%     Which reasonably approximates the analytical global minimum (note
%     that due to randomness, your results will likely not be exactly the
%     same).

%  Reference:
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.

%   joachim.vandekerckhove@psy.kuleuven.be
%   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $


InitS= 4;  % Initial Step Size.  Step ranges from -1 to +1  or +-5 %4

StepDecayFactor= 0.98;
TempDecayFactor= 0.98;
InitT= 1;

def = struct(...
        'CoolSched',@(T) (TempDecayFactor*T),...  
        'BasinFitness', @(x) sum((bias*round(abs(x)))) , ...  % relative frequency:  sum(- 15/d* (cos(2*pi*x)-1)) % Local optima fitness:  sum((3*round(abs(x))))
        'Generator',@(p,q) (p + q), ...
        'StepSize',@(x,tc) ((rand(1, length(x)).*(InitS*2))-InitS).*(StepDecayFactor^tc),...  %%rand uniform num in any direction add or subtracy (rand num -1 and +1 then add to x)   % x+((rand(1, length(x)).*2)-1)????   % (((rand(1, length(x)).*2.5)-1.25))
        'MaxStepSize', @(x,tc) ones(1, length(x)).*(max(abs(([0,1].*(InitS*2))-InitS))).*(StepDecayFactor^tc),...    %ones becuase rand's max value is 1
        'InitTemp',InitT,...
        'MaxConsRej',20000,...  % 1000% of MaxTries%2400
        'MaxSuccess',60,...   %  25% of MaxTries
        'MaxTries',240,...    % controls constant nature of step size plot
        'StopTemp',1e-8,...
        'StopVal',-Inf,...
        'StopStepVal', Inf,...   % Step at which the search is stopped
        'Verbosity',1);

% Check input
if ~nargin %user wants default options, give it and stop
    minimum = def;
    return
elseif nargin<2, %user gave only objective function, throw error
    error('MATLAB:anneal:noParent','You need to input a first guess.');
elseif nargin<3, %user gave no options structure, use default
    options=def;
else %user gave all input, check if options structure is complete
    if ~isstruct(options)
        error('MATLAB:anneal:badOptions',...
            'Input argument ''options'' is not a structure')
    end
    fs = {'CoolSched','BasinFitness','Generator','StepSize', 'MaxStepSize', 'InitTemp','MaxConsRej',...
        'MaxSuccess','MaxTries','StopTemp','StopVal', 'StopStepVal', 'Verbosity'};
    for nm=1:length(fs)
        if ~isfield(options,fs{nm}), options.(fs{nm}) = def.(fs{nm}); end
    end
end

% main settings
newsol = options.Generator;      % neighborhood space function
basinFitness= options.BasinFitness;
Tinit = options.InitTemp;        % initial temp
minT = options.StopTemp;         % stopping temp
cool = options.CoolSched;        % annealing schedule
Sinit= options.StepSize;         %Initial Step size
max_stepSize= options.MaxStepSize;       % The maximum step size for each iteration
minF = options.StopVal;
max_consec_rejections = options.MaxConsRej;
max_try = options.MaxTries;
max_success = options.MaxSuccess;
stop_StepValue= options.StopStepVal;
report = options.Verbosity;
k = 1;                           % boltzmann constant

% counters etc

% Tracking All new points (rejections)
all_counter=0;  

fval_all= ones(20000, 1)* (-1);


% Tracking Accepted points
acc_counter=0;  %%accepted counter
fval_accepted= ones(20000, 1)* (-1);
accepted_itertaions= ones(20000, 1)* (-1);



%%% f(x) tracker
counter=1;  %new  % f(x) counter
x_tracker= ones(20000,d)*50;    % initilizing 
x_tracker(counter,:)= parent;    


itry = 0;
success = 0;
finished = 0;
consec = 0;
T = Tinit;
initenergy = loss(parent);
oldenergy = initenergy;



fval_tracker= ones(20000, 1)* (-1);    % initilizing 
fval_tracker(counter)= oldenergy;      

total = 0;
if report==2, fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy); end


% To count rejections in first 100 ietrations
rej_counter=0;  
Rej_rate= ones(20000, 1)* (-1);

%Counters to track temp chnages
temp_counter= 0; % Array index for Temp_change_iterations
                 


% Counter to track backward step acceptance rate at each step or iteration
backward_acc_counter=0;  %%Counts how many backward steps were accepted
last_BackAccRate= ones(20000, 1)* (-1); % stores backward step accpetance rate for each last 100 iterations

                 
                 
%%% Step Size Analysis
                 
csr=0;  % counter for step size red which are improving & accepeted point
sr= ones(20000, 1)* (30);
isr= ones(20000, 1)* -1;


csb=0;  % counter for step size red which are improving & accepeted point
sb= ones(20000, 1)* (30);
isb= ones(20000, 1)* -1;

csblue=0;
sblue= ones(20000, 1)* (30);
isblue= ones(20000, 1)* -1;

Max_S= ones(20000, 1)* -1;



%%Relative Frequency of accepted values only
rf_acc = ones(20000, 1)* (-100);   % same as n_trial_best in random walk exp % initilizing relative frequency= newergy(of x)- sum((bias*round(abs(x))))


%accepted backward steps: 
%count of accepted steps when basinFitness(newenergy) > basinFitness(oldenergy) 
back_acc_counter=0;



                 
while (InitS*StepDecayFactor^temp_counter)>0.1 
    itry = itry+1; % just an iteration counter
    current = parent; 
    S= Sinit(current, temp_counter);  %new  % inital and all other Step Size
    maxS= max_stepSize(current, temp_counter);
    Max_S(counter)=  sqrt(sum(maxS.^2));   % sqrt(sum((maxS.^2))); % caldulates length of step
    
    counter= counter+1;   % f(x) counter
    
    % % Stop / decrement T criteria  %Decrement Step size S
    if itry >= max_try || success >= max_success;
        if T < minT || consec >= max_consec_rejections;
            finished = 1;
            total = total + itry;
            break;
        else
           
            T = cool(T);  % decrease T according to cooling schedule
            temp_counter= temp_counter+1;
            
            
            if report==2, % output
                fprintf(1,'  T = %7.5f, loss = %10.5f\n',T,oldenergy);
            end
            total = total + itry;
            itry = 1;
            success = 1;
        end
    end
    
    newparam = newsol(current, S);  % Adds current to the step size  %Generator= current+S
    newenergy = loss(newparam);
    
    all_counter= all_counter+1; 
    fval_all(all_counter)= newenergy- basinFitness(newparam);  % newenergy-relativeFrequency(newparam);
    
    
    %{
    %% This is only executed if you have a stop value or minF (our case is -Inf)
    
    if (newenergy < minF),  
        parent = newparam; 
        x_tracker(counter, :)= parent; 
        oldenergy = newenergy;
        fval_tracker(counter)=oldenergy; 
        
        acc_counter= acc_counter+1; 
        x_accepted(acc_counter, :)= parent; 
        fval_accepted(acc_counter)=  oldenergy- relativeFrequency(parent) ; %new % oldenergy- relativeFrequency(parent)
        accepted_itertaions(acc_counter)= all_counter; 
        
        
        break
    end
    %}
    
    parent_basinFitness= basinFitness(parent); 
   
    
    if (oldenergy-newenergy > 1e-6)
        
        parent = newparam;
        x_tracker(counter,:)= parent; 
        oldenergy = newenergy;
        fval_tracker(counter)= oldenergy; 
        success = success+1;
        consec = 0;     
       
        acc_counter= acc_counter+1; 
   
        fval_accepted(acc_counter)=  oldenergy- basinFitness(parent); %new % oldenergy- relativeFrequency(parent)
        accepted_itertaions(acc_counter)= all_counter;  
        
        csr= csr+1;    % No. of steps for improving steps that are accepted
        sr(csr)= sqrt(sum((S.^2)));  % Red Step Size which indicates improving steps that are accepted
        isr(csr)= all_counter;   %step number at which they occur
        
        rf_acc(acc_counter) =  oldenergy- basinFitness(parent); %relative frequency of accepted (improving and disimproving) values only % Excludes the first relative frequency but includes the last one (note last not not used in iteration!)
        
        
        % Counts the accepetd step that are deceptive or backward steps (as in RW)
        if (basinFitness(newparam)> parent_basinFitness)
            back_acc_counter= back_acc_counter+1;  
        end
        
    else
        if (rand < exp( (oldenergy-newenergy)/(k*T) ));    % deltaE= (oldenergy-newenergy) find a T for wjich exp will be 0.5
            parent = newparam;
            x_tracker(counter, :)= parent; 
            oldenergy = newenergy;
            fval_tracker(counter)=oldenergy; 
            success = success+1;
           
            acc_counter= acc_counter+1; 
                    %# x_accepted(acc_counter, :)= parent; 
            fval_accepted(acc_counter)=  oldenergy-basinFitness(parent) ; % oldenergy- relativeFrequency(parent)
            accepted_itertaions(acc_counter)= all_counter; 
            
            backward_acc_counter= backward_acc_counter+1;  % counts backward steps that are accepted
            
            csb= csb+1;    % No. of steps for disimproving steps that are accepted
            sb(csb)= sqrt(sum((S.^2)));  % Red Step Size which indicates disimproving steps that are accepted
            isb(csb)= all_counter;  %step number at which they occur
            
            rf_acc(acc_counter) =  oldenergy- basinFitness(parent); %relative frequency of accepted (improving and disimproving) values only % Excludes the first relative frequency but includes the last one (note last not not used in iteration!)
            
            
            % Counts the accepetd step that are deceptive or backward steps (as in RW)
            if (basinFitness(newparam)> parent_basinFitness)
                back_acc_counter= back_acc_counter+1;  
            end
            
            
        else
            consec = consec+1;
           
            x_tracker(counter, :)= parent; 
            fval_tracker(counter)=oldenergy; 
            
            % Counting no. of rejections in the first 100 iterations only
            if (counter<=100)          
                rej_counter=rej_counter+1;  
            end                    
            
            csblue= csblue+1;    % No. of steps for disimproving steps that are rejected
            sblue(csblue)= sqrt(sum((S.^2)));  % Red Step Size which indicates disimproving steps that are accepted
            isblue(csblue)= all_counter; %step number at which they occur
        end
    end
    
    last_BackAccRate(all_counter)  = (backward_acc_counter/all_counter)*100;
end

minimum = parent;
fval = oldenergy; %fval is the value of the LOSS function evaluated at MINIMUM (with minimum as input variable).

xTracker= x_tracker;
fvalTracker= fval_tracker;

rejs= rej_counter;  % Rejections for first 100 iterations only

All_energies= fval_all; %Relative frequency of all energies of the all positions


acc_energies= fval_accepted;  % Relative frequency of accepted positions
acc_it= accepted_itertaions; 

%Tracking temp chnage
last_BAccR=last_BackAccRate;

% Tacking Step Size
sR= [isr, sr];
sBlk= [isb, sb];
sBlue= [isblue, sblue];


% Tracking Relative Frequency of accepted values
rf_accepted=rf_acc;

%Counting accepted backward (or deceptive) steps
back_acc_count=back_acc_counter;

if report;
    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', total);
    fprintf(1, '  Total final loss:        \t%g\n', fval);
end

