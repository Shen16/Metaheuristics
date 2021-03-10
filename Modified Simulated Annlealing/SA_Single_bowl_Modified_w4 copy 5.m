%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   anneal.m  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Usage:
%     [MINIMUM,FVAL] = ANNEAL(LOSS,NEWSOL,[OPTIONS]);
%          MINIMUM is the solution which generated the smallest encountered
%          value when input into LOSS.
%          FVAL is the value of the LOSS function evaluated at MINIMUM.



%%%%%%%
% Rastrigin function handle with anneal
d=1;

basin_height=600;
bias=3;
initPos= 500.25;

rastrigin =  @(x)sum( (bias*round(abs(x))) - (basin_height/2)/d* (cos(2*pi*x)-1)); %x is a row vector of length d


%Intial position parent
%{
a = -initPos;
b = initPos;
parent = (b-a).*rand(1,d) + a; % Generates a random intial position of d dimension, within +5.12 and -5.12 range
%}

parent = ones(1,d)* initPos; % Starts from 500.25 or 5.25

[x, f, rej, xT, fT, Af, Sf, Si ,BAccR ,sR ,sBlk ,sBlue, maxS, minS] = anneal(rastrigin, parent, struct('Verbosity', 1), d, bias );  % x= final position, f= final energy, rej= rejections in first 100 iterations only, xT= position tracker, fT= energy tracker, Ax= All new positions and Af= All new solution energies, Ai= All solution iterations, Sx= Position accepted, Sf= Energies of accepted positions, Si= A vector containing step number at which acceptance occured, temp_it= iterations at which temp chnages





%%%%%%%%%%%%%%%%%%%%%
% 1,2, or 3D Plot
figure;
if size(xT,2)==1
    % 1D plot
    x= xT(xT~=50,1);
    plot(1:length(x),x)
    ylabel('Position of X')
    xlabel('Number of Iterations')
    title('1D Simulated Annealing')
    
elseif size(xT,2)==2
    % 2D plot
    x= xT(xT(:,1)~=50,1);
    y= xT(xT(:,2)~=50,2);
    plot(x,y,'.','MarkerEdgeColor','r', 'LineStyle', '-')
    
    hold on;
    grid on %grid lines
    xlabel('X','fontSize',12);
    ylabel('Y','fontSize',12);
    title('2D Simulated Annealing','fontsize',14)
    plot(0,0, 'Marker', '*', 'Color', 'black', 'MarkerSize' , 10)  %Plotting origin
    
elseif size(xT,2)==3
    % 3D Plot of x
    x= xT(xT(:,1)~=50,1);
    y= xT(xT(:,2)~=50,2);
    z= xT(xT(:,3)~=50,3);
    
    plot3(x,y,z,'.','MarkerEdgeColor','r', 'LineStyle', '-')
    
    hold on;
    grid on %grid lines
    view(3) %view 3D
    xlabel('X','fontSize',12);
    ylabel('Y','fontSize',12);
    zlabel('Z','fontSize',12);
    title('3D Simulated Annealing','fontsize',14)
    plot3(0,0,0, 'Marker', '*', 'Color', 'green', 'MarkerSize' , 10)
else
    disp("Cannot create a plot for this dimension!");
end




%%%%%%%%%%%%%%%%%%%%%
% Plot of f(x) against number of iterations
figure;
f= fT(fT~=-1);
iterations= 1:length(f);
plot(iterations, f, '.','MarkerEdgeColor','r', 'LineStyle', '-')
xlabel('Step Count','fontSize',12);
ylabel('f(x)','fontSize',12);
t1= sprintf('Change of f(x) in %d dimensions', d);
title(t1)



%%%%%%%%%%%%%%%%%%%%%
% Distance from origin for each x position
figure;
d= size(xT,2);
for i=1:d
    v(:,i)= xT(xT(:,i)~=50,i);
end
v;

origin=zeros(1,d);
rows= size(v,1);
dist_frm_origin= zeros(rows, 1);
for j=1:rows
    dist_frm_origin(j)= sqrt(sum((v(j,:)-origin).^2));
end

%plot against number of iterations
iterations= 1:length(dist_frm_origin);
plot(iterations, dist_frm_origin, '.','MarkerEdgeColor','r', 'LineStyle', '-')
xlabel('Step Count','fontSize',12);
ylabel('Distance of position x from Origin','fontSize',12);
t2= sprintf('Distance from Origin in %d dimensions', d);
title(t2)



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



%%% A plot of the percentage of accepted backward steps at each step
BackAccRate= BAccR(BAccR~=-1);

StepCount= 1:length(BackAccRate);

figure;
plot(StepCount, BackAccRate);

length(StepCount);


xlabel('Step Count','fontSize',12);
ylabel('Back Acceptance Rate','fontSize',12);
t1= sprintf('Back Acceptance Rate for each step in %d dimensions', d);
title(t1);




%%% Actual Step Size plot
figure;
hold on;

Stepit_Blue= sBlue(sBlue(:,1)~=-1,1);
StepSize_Blue= sBlue(sBlue(:, 2)~=30,2);
plot(Stepit_Blue, StepSize_Blue, '.b');

Stepit_Red= sR(sR(:,1)~=-1,1);
StepSize_Red= sR(sR(:, 2)~=30,2);
plot(Stepit_Red, StepSize_Red, 'or', 'Markersize', 6);

Stepit_Black= sBlk(sBlk(:,1)~=-1,1);
StepSize_Black= sBlk(sBlk(:, 2)~=30,2);
plot(Stepit_Black, StepSize_Black, '.k', 'Markersize', 8);

%{
Max_StepSize= maxS(maxS~=-1);
f= fT(fT~=-1);
all_iterations= length(f);
plot(1:all_iterations, Max_StepSize, 'lineStyle', '-', 'Color', 'cyan');

Min_StepSize= minS(minS~=-1);
plot(1:all_iterations, Min_StepSize, 'lineStyle', '-', 'Color', 'cyan');
%}

xlabel('Step Count','fontSize',12);
ylabel('Step Size','fontSize',12);
t1= sprintf('Step Size Analysis in %d dimensions', d);
title(t1);


hold off;


