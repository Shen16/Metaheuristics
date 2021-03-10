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
d=1;
n=30;   % Number of simulations

basin_height=600;
bias=3;
initPos= 500.25;

x= zeros(n,d);  % To record the final psoition x (of d dimension) after each simulation
f= zeros(n,1);  % To record the final f(x) after each simulation
rej= zeros(n,1);
dist= zeros(n,1);


n_mean_rf_accs= zeros(n,1); %NOT NEEDED anymore % To record the mean values of rf_acc values for each n iterations

back_acc_count= zeros(n,1); %Storing counts of accepted backward (or deceptive) steps for n iterations



for i=1:n
    
    rastrigin = @(x) sum( (bias*round(abs(x))) - (basin_height/2)/d* (cos(2*pi*x)-1)); % x is a row vector of length d  %% make cos part 0 to 30 again (decides depth)  % global (0,0)
    
    %Intial position parent
    %{
    a = -initPos;
    b = initPos;
    parent = (b-a).*rand(1,d) + a; % Generates a random intial position of d dimension, within +5.12 and -5.12 range
    %}
    
    parent= ones(1,d)* initPos; % Starts from 500.25 (with InitS=5) or 5.25 (with InitS=1)
    
    
    
    [p, q, r, xT, fT, Af, Sf, Si, BAccR, sR, sBlk, sBlue, maxS, rf_acc, back_acc] = anneal(rastrigin, parent, struct('Verbosity', 0, 'StopTemp',1e-8, 'StopVal',-Inf), d ,bias);
    x(i,:)= p;
    f(i)=q;  %Last fitness value of each iteration
    rej(i)=r;
    
    % Initial Step Size
    dist(i)= sqrt(sum((xT(2,:)-xT(1,:)).^2));
    
    
    %Relative Frequency for accepted points only
    rf_accepted= rf_acc(rf_acc~=-100);
    n_mean_rf_accs(i)= mean(rf_accepted); % mean of rf_acc values for each n iterations
    
    %Counts of accepted backward (or deceptive) steps
    back_acc_count(i)= back_acc;
    
    
    
end


initial_step_size = mean(dist);

time= toc; %end timer


%{

name= sprintf('%dD_SA_new', d);
filename= name + ".mat";
save(filename, 'f','d', 'x','rej','time', 'initial_step_size','p', 'q', 'r', 'xT', 'fT', 'Af', 'Sf', 'Si', 'BAccR', 'sR', 'sBlk', 'sBlue', 'maxS', 'n_mean_rf_accs', 'rf_accepted', 'back_acc_count')

%}



%%Load results to plot
load('30D_SA_new.mat');
%n=30

%To check intial plot
%figure;
%histogram(f, 1000)

%%%%%%%%%%%%%%%%%%%%%
% Histogram and mean of f(x)


x_min= 0; % x-axis range minimum %-10
x_max= 36000; %26000
y_min=0;
y_max= n;  %14000


edges = linspace(x_min, x_max, 200); % Create 20 bins
figure;
histogram(f,'BinEdges',edges)  %, 'Binwidth', 0.1



xlim([x_min, x_max]);
ylim([y_min, y_max])


%set(gca,'XTick',[], 'YTick', [])



%%Remove tick labels for the X and Y axes
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

%%To remove tick marks on the y-axis
%tickMarks = {'YTick',[]};
%set(gca,tickMarks{:});

%ticks from -100 to 100 at intervals of 100
xticks(x_min:2000:x_max )
yticks(y_min:10:y_max)




%xlabel('f(x) values','fontSize',12);
%ylabel('Relative Frequency of f(x) value','fontSize',12);
%t1= sprintf('Histogram of f(x) in %d dimensions', d);
%title(t1)

mf=mean(f);  %mean of f(x)
sdf= std(f); %sd of f(x)

%mean of backward steps or rejections
mr= mean(rej);


%Relative Frequency for accepted points only --> mean of n (=100) cycles
mrf=mean(n_mean_rf_accs);


%Counts of accepted backward (or deceptive) steps
mean_backAcc= mean(back_acc_count);

T= table(d,mf, sdf, time, mr, initial_step_size, mrf, mean_backAcc);
T.Properties.VariableNames = { 'Dimension' , 'Mean of f(x)' , 'Sd of f(x)', 'Time Elapsed', 'Mean of Rejections', 'Mean Initial Step Size' , 'Mean Relative Frequency', 'Mean Backward steps'}


%{

%%%%%%%%%%%%%%%%%%%%%
% Histogram of Distance from origin for each x position
%{
figure;
x_min= 0; % x-axis range minimum %-10
x_max= 7000; %600
y_min=0;
y_max= 30;  %14000
%}

d= size(x,2);

origin=zeros(1,d);
rows= size(x,1);
dist_frm_origin= zeros(rows, 1);
for j=1:rows
    dist_frm_origin(j)= sqrt(sum((x(j,:)-origin).^2));
end

mo= mean(dist_frm_origin);  %mean of distnace from origin
sdo= std(dist_frm_origin);  %standard deviation of distnace from origin

%{
histogram(dist_frm_origin, 'BinEdges',edges)   %'BinWidth', 0.1, 'EdgeColor', 'red'
edges = linspace(x_min, x_max, 10000); % Create 20 bins



xlim([x_min, x_max]);
ylim([y_min, y_max])
%set(gca,'XTick',[], 'YTick', [])


%%Remove tick labels for the X and Y axes
set(gca,'Yticklabel',[])
%set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

%%To remove tick marks on the y-axis
tickMarks = {'YTick',[]};
set(gca,tickMarks{:});

ticks= 0:1000:x_max  %ticks from -100 to 100 at intervals of 100
xticks(ticks)

%}

%To check intial plot
figure;
histogram(dist_frm_origin, 1000)

xlabel('Distance of x from Origin','fontSize',12);
ylabel('Relative Frequency','fontSize',12);
t2= sprintf('Distance from Origin in %d dimensions', d);
title(t2)



%mean of backward steps or rejections
mr= mean(rej);

T= table(d,mf, sdf, time, mr, initial_step_size);
T.Properties.VariableNames = { 'Dimension' , 'Mean of f(x)' , 'Sd of f(x)', 'Time Elapsed', 'Mean of Rejections', 'Mean Initial Step Size'}


%{
diary off

type 1D_SA_Rastr
type 2D_SA_Rastr
type 3D_SA_Rastr
type 5D_SA_Rastr
type 10D_SA_Rastr
type 30D_SA_Rastr

delete('30D_SA_Rastr')

%}


%}

%load('30D_SA.mat');
mf=mean(f)  %mean of f(x)

d= [1,2,3,5,10,20,30]
mean_list= [0.5100, 2.5642e+03, 4.3498e+03, 7.4349e+03, 1.4799e+04, 2.6726e+04, 2.9463e+04]

figure;
plot(d, mean_list)
title("Mean Fitness")
legend("Normal SA", "Location", "southeast")

x_min= 0; % x-axis range minimum %-10
x_max= 30; %26000
y_min=0;
y_max= 30000;  %14000


xlim([x_min, x_max]);
ylim([y_min, y_max])


%%%%%%%%%%%%%%%%%%%%%
% Plot Relative Frequency of f(x). Accepted (red) and rejected (blue) against number of iterations
figure;
All_energies= Af(Af~=-1);
relative_frequency_all= All_energies;
iterations= 1:length(All_energies);
plot(iterations, relative_frequency_all, '.','MarkerEdgeColor','b', 'LineStyle', 'none', 'MarkerSize', 0.1)
xlabel('Step Count','fontSize',12);
ylabel('Relative Frequency of f(x)','fontSize',12);
t1= sprintf('Relative frequency of f(x) in %d dimensions', d);
title(t1)


hold on;
Succesful_energy= Sf(Sf~=-1);
relative_frequency_s= Succesful_energy;
succesful_iterations= Si(Si~=-1);
plot(succesful_iterations, relative_frequency_s , '.','MarkerEdgeColor','r', 'LineStyle', 'none', 'MarkerSize', 10);

hold off;

All_energies
mean(All_energies)



%{
%Relative Frequency for best points only--> Wrong since included all points, not just the accepetd ones

load('30D_SA.mat');



basin_height=600;
bias=3;

rastrigin = @(x) sum( (bias*round(abs(x))) - (basin_height/2)/d* (cos(2*pi*x)-1));
basinFitness= @(x) sum((bias*round(abs(x)))) ;


d= size(xT,2);
for i=1:d
    v(:,i)= xT(xT(:,i)~=50,i);
end
v;

rows= size(v,1);
relativeFrequency= zeros(rows, 1);
for j=1:rows
    relativeFrequency(j)= rastrigin(xT(j))- basinFitness(xT(j)); 
end


figure;

%plot against number of iterations
iterations= 1:length(relativeFrequency);
plot(iterations, relativeFrequency, '.','MarkerEdgeColor','r', 'LineStyle', '-')
xlabel('Step Count','fontSize',12);
ylabel('Relative Frequency','fontSize',12);
t2= sprintf('Relative Frequency of Best position in %d dimensions', d);
title(t2)

mean(relativeFrequency(1:100))


Dimension = [1;2;3;5;10;20;30];
Mean = [10.4429;14.6070;5.8212;32.9155;1.2476;14.6600;8.1315];


T = table(Dimension, Mean)

%}

%Relative Frequency for accepted points only --> mean of n (=100) cycles
mean(n_mean_rf_accs)


Dimension = [1;2;3;5;10;20;30];
Mean_rf = [1.3306;6.9929;22.095;55.623;107.59;195.18;244.53];


SA_normal = table(Dimension, Mean_rf)

% Count of Backward accepted steps
Dimension = [1;2;3;5;10;20;30];
Mean_back_acc = [25.067 ;8.0667 ; 5.1333 ;4.1667; 4.9 ; 22.167; 86.833];


SA_normal = table(Dimension, Mean_back_acc)



