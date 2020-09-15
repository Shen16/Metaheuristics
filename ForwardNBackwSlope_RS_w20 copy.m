% Forward & Backward with slope Random Search (FBS_RS)
% Minimization Problem
% n=10000, m=100000


%diary 1D_Uni_FBS

%%%%%%%%%%%%%%%%%%%%%%%%
%%% OneD_Uniform Random Walk
%%%%%%%%%%%%%%%%%%%%%%%%

tic; %start timer

d=1;
basin_height= 600;

m=100000;
finalPos= zeros(1,m); %for m iterations below 
finalbest= zeros(1,m); % to track final best values

totalSuccessful= zeros(1,m); 
totalDeceptive= zeros(1,m); 


for i=1:m 
    
    n=10000;  
    position=0;
    
    best= basin_height;  %new % previously best=30, now best=600
    successful=0;
    deceptive=0;
    
    
    for j=1:n
        r1= rand* basin_height;   %new % prev rand*30
        r2= rand* basin_height;   %new

        if (r1-3)< best && (r1-3)< (r2+3)   %Making r1 the better attarction basin (minimization problem) and r2 the worse
            position= position+1;
            successful= successful+1;
            best= r1;
        elseif (r2+3)< best && (r2+3)< (r1-3)
            position= position-1;
            deceptive = deceptive+1;
            best=r2;
        else
            position=position;

        end
    end
    finalPos(i)= position;
    totalSuccessful(i)= successful;   
    totalDeceptive(i)= deceptive;   
    finalbest(i)= best;
end
time =toc; %end timer

%Poisson Distribution
histogram(finalPos)  %histogram(result), result== final pos of tracker
ylabel('Relative Freqeuncy of Position')
xlabel('Position')
title('1D Uniform Random Search Forward & Back with Slope')
dim = [.2 .5 .3 .3];
str = [strcat("m (outerloop)=",num2str(m));strcat( "n (innerloop)=",num2str(n))];
annotation('textbox',dim,'String',str,'FitBoxToText','on')

%{
data = finalpos %0.1 * rand(1, 1000); % Create sample data. 
edges = linspace(0, 1, 21); % Create 20 bins.
% Plot the histogram.
histogram(data, 'BinEdges',edges);
% Fancy up the graph.
grid on;
xlim([0, 1]);
xlabel('Data Value', 'FontSize', 14);
ylabel('Bin Count', 'FontSize', 14);
title('Histogram of Data', 'FontSize', 14);

%}

%Statistics
format longG
Mean= mean(finalPos);
variance= var(finalPos);
sd= std(finalPos);
maxFinalPos= max(finalPos);
minFinalPos= min(finalPos);

successful_mean= mean(totalSuccessful);  
deceptive_mean= mean(totalDeceptive); 

best_mean= mean(finalbest); 

T1= table(d,Mean,variance, sd, maxFinalPos, minFinalPos, time, successful_mean, deceptive_mean, best_mean);
T1.Properties.VariableNames = { 'Dimension' , 'Mean' , 'Variance', 'Standard Dev', 'Max Pos', 'Min Pos', 'Time', 'Mean Successful', 'Mean Deceptive', 'Mean Best'}

%diary off

%type 1D_Uni_FBS  % fig: 1D_Uni_FBS (for n=10,000)


%log(ms+md)

%delete("1D_Uni_FBS")







%diary 30D_Cos_FBS
%%%%%%%%%%%%%%%%%%%%%%%%
%%% All_D_Cos
%%%%%%%%%%%%%%%%%%%%%%%%

tic; %start timer
d=30;
basin_height= 600;

m=100000;
finalPos= zeros(1,m); %for m iterations below 
finalbest= zeros(1,m); % to track final best values

totalSuccessful= zeros(1,m);  
totalDeceptive= zeros(1,m);  



for i=1:m 
    
    n= 10000;
    position= 0; 
    
    best=basin_height;
    successful=0;
    deceptive=0;
    
    rand_1= zeros(1,n); %vector. A rand number for each simulation
    rand_2= zeros(1,n);

    for j=1:n

        % Initialize first set of random numbers : rand1

        value1=0;
         for k=1:d
             value1= value1 + (((basin_height/2)/d)*cos(2*pi*rand));
         end
        rand_1(i)=value1+ (basin_height/2);


        % Initialize second set of random numbers : rand2

        value2=0;
         for k=1:d
             value2= value2 + (((basin_height/2)/d)*cos(2*pi*rand));  
         end

        rand_2(i)=value2+ (basin_height/2);
        
        %Move forward if first random number (from better attraction basin)
        %is less than best and less than second random number (from worse attraction basin)
        if (rand_1(i)-3)< best && (rand_1(i)-3)< (rand_2(i) +3)
            position=position +1;
            successful= successful+1;
            best= rand_1(i);
        elseif (rand_2(i) +3) <best && (rand_2(i) +3) < (rand_1(i)-3)
            position= position-1;
            deceptive = deceptive + 1;
            best = rand_2(i);
        else
            position=position;
        end
    end
    finalPos(i)= position;
    totalSuccessful(i)= successful;
    totalDeceptive(i)= deceptive;
    finalbest(i)= best;
    
end
time= toc; %end timer
    
%Poisson Distribution
figure ;
histogram(finalPos)
xlabel('Position') ;
ylabel('Relative Frequency of Position') 

t= sprintf('Random Search Forward & Back with Slope in %d dimensions', d);
title(t)

dim = [.2 .5 .3 .3];
str = [strcat("m (outerloop)=",num2str(m));strcat( "n (innerloop)=",num2str(n))];
annotation('textbox',dim,'String',str,'FitBoxToText','on')

    
%Statistics
Mean= mean(finalPos);
variance= var(finalPos);
sd= std(finalPos);
maxFinalPos= max(finalPos);
minFinalPos= min(finalPos);

successful_mean= mean(totalSuccessful);  
deceptive_mean= mean(totalDeceptive);  

best_mean= mean(finalbest); 

T2= table(d,Mean,variance, sd, maxFinalPos, minFinalPos, time, successful_mean, deceptive_mean, best_mean);
T2.Properties.VariableNames = { 'Dimension' , 'Mean' , 'Variance', 'Standard Dev', 'Max Pos', 'Min Pos', 'Time', 'Mean Successful', 'Mean Deceptive', 'Mean Best'}

%{
diary off

type 1D_Cos_FBS % fig: 1D_Cos_FBS
type 2D_Cos_FBS % fig: 2D_Cos_FBS
type 3D_Cos_FBS % fig: 3D_Cos_FBS
type 5D_Cos_FBS % fig: 5D_Cos_FBS
type 10D_Cos_FBS % fig: 10D_Cos_FBS
type 30D_Cos_FBS % fig: 30D_Cos_FBS



%Growth of function for n=10000

x= [1,2,3,5,10,30]  % for cosine dist
y= [185.11, 264.25, 412.23, 638.44, 1266.05, 3566.75]
bigOh= scatter(x,y)
lsline
title("Running Time Complexity of Forward & Back with Slope Random Search")
xlabel("Number of Dimension")
ylabel("Running Time")




%Mean Position or solution produced as dimensions increase for n=10000
x= [1,2,3,5,10,30]  % for cosine dist
y= [2347.29, 1142.01, 816.13, 889.79, 1999.56, 6727.61] % Mean of Successful steps
bigOh= plot(x,y)
title("Successful Step Rate of Forward & Back with Slope Random Search")
xlabel("Number of Dimension")
ylabel("Mean Sucessful Steps")   % the lower the better since minimization problem

%}


%{
%Same plot as above
%Rate of Success as dimensions increase for n=10000
figure;
x= [1,2,3,5,10,30]  % for cosine dist
y= [2415.85, 1210.52, 876.99, 955.15, 2094.46, 6776.57] % Mean of Successful steps
bigOh= plot(x,y)
title("Running Time Complexity of Forward & Back with Slope Random Search")
xlabel("Number of Dimension")
ylabel("Successful Steps")

%}


%delete("1D_Cos_FBS")

