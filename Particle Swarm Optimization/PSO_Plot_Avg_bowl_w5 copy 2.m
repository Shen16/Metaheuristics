%%%%%% PSO Plots only

load("30D_PSO_Rastr.mat")

%%%%%%%%%%%%%%%%%%%%%
% Histogram and mean of f(x)

%{
% bin specs.
nbins = 50;
bound = 150;
bins = linspace(0,bound,nbins);
%}

x_min= 0; % x-axis range minimum %-10
x_max= 220; %600
y_min=0;
y_max= 100;  %14000

edges = linspace(x_min, x_max, 1000); % Create 20 bins
figure;
histogram(fval, 'BinEdges',edges)  %, 'Binwidth', 0.1

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


mf=mean(fval);  %mean of f(x)
sdf= std(fval); %sd of f(x)


%xlabel('NewCost values','fontSize',12);
%ylabel('Relative Frequency of f(x) or NewCost value','fontSize',12);
%t1= sprintf('Histogram of f(x) in %d dimensions', d);
%title(t1)

%Counts of accepted backward (or deceptive) steps
mean_backAcc= mean(back_acc_count)
format shortG


T= table(d,mf, sdf, time,mean_backAcc );
T.Properties.VariableNames = { 'Dimension' , 'Mean of f(x)' , 'Sd of f(x)', 'Time Elapsed', 'Mean Backward steps'}




% Count of Backward accepted steps
Dimension = [1;2;3;5;10;20;30];
Mean_back_acc = [985.33; 2464.6; 3655; 5224.2; 6902.2; 8014.9 ;8209.8];


PSO_normal = table(Dimension, Mean_back_acc)



%{

%%%%%%%%%%%%%%%%%%%%%
% Histogram of Distance from origin for each x position
figure;

origin=zeros(1,d);
dist_frm_origin= zeros(n, 1); % n=100 % No . of simulation

for j=1:n
    dist_frm_origin(j)= sqrt(sum((gbest_x(j,:)-origin).^2));
end

histogram(dist_frm_origin, 1000)   %'BinWidth', 0.1, 'EdgeColor', 'red'
mo= mean(dist_frm_origin);  %mean of distnace from origin
sdo= std(dist_frm_origin);  %standard deviation of distnace from origin

xlabel('Mean Distance from origin of the particle at last step that has the minimum fitness value','fontSize',12);
ylabel('Relative Frequency','fontSize',12);
t2= sprintf('Distance from Origin in %d dimensions', d);
title(t2)




T= table(d,mf, sdf, time);
T.Properties.VariableNames = { 'Dimension' , 'Mean of f(x)' , 'Sd of f(x)', 'Time Elapsed'}

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


