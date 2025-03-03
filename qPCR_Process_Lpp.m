%% loading and processing data from qPCR
% My goal in this code is to load the raw qPCR data from my experiments,
% specifying the cells of the excel file to use for the data and for the
% data structure name using 'readcell'. 

% Once the data is loaded, each individual dataset is run through a
% function that calculates the proper baseline value from the slopes in the
% top and bottom half of log(data) before the plateau based on the
% algorithm in the paper:
% "Amplification efficiency: linking baseline and bias in the analysis of quantitative PCR data"

% Once that is finished, the efficiency of each sample is calculated by
% finding the slope of the line before the plateau. The plateau begins when
% the the second derivative of the raw data is at a minimum.

%%
clear
%%
% Exp 1 is 2022.08.22, and has A-F 1-8,10
% Exp 2 is 2022.12.21, and has A-H 1-12
% Exp 3 is 2023.11.14, and has A-F 1-7
% Exp 4 is 2023.12.5, and has A-G 1-9

FILENAME = "/Users/bpweaver/Desktop/Research/LppProject/Raw_qPCR.xlsx"; 
% Check if this has been moved to another folder and update this if
% necessary.

%% Variables for naming parts of data structures
% Letters and numbers here correspond to wells in the qPCR machine,
% and individual data sets in each experiment.
letters1 = 'A':'F';
letters2 = 'A':'H';
letters3 = 'A':'F';
letters4 = 'A':'G';

list1 = [1:8,10];
list2 = 1:12;
list3 = 1:7;
list4 = 1:9;

% Loading data into structures
qPCR_data.exp1 = qPCR_load(FILENAME,"2022.08.22",letters1,list1);
qPCR_data.exp2 = qPCR_load(FILENAME,"2022.12.21",letters2,list2);
qPCR_data.exp3 = qPCR_load(FILENAME,"2023.11.14",letters3,list3);
qPCR_data.exp4 = qPCR_load(FILENAME,"2023.12.5",letters4,list4);

%% Calculating Baseline and subtracting from the data
BC_qPCR.exp1 = Baseline_calc(letters1, list1, qPCR_data.exp1);
BC_qPCR.exp2 = Baseline_calc(letters2, list2, qPCR_data.exp2);
BC_qPCR.exp3 = Baseline_calc(letters3, list3, qPCR_data.exp3);
BC_qPCR.exp4 = Baseline_calc(letters4, list4, qPCR_data.exp4);

% Calculating individual efficiencies
qPCR_eff.exp1 = Efficiencies(letters1,list1,BC_qPCR.exp1);
qPCR_eff.exp2 = Efficiencies(letters2,list2,BC_qPCR.exp2);
qPCR_eff.exp3 = Efficiencies(letters3,list3,BC_qPCR.exp3);
qPCR_eff.exp4 = Efficiencies(letters4,list4,BC_qPCR.exp4);

%% Grouping data from each experiment based on conditions
cond1 = ["WT","D","0ng","10ng","50ng","100ng"];
cond2 = ["WT","D","0ng","10ng","50ng","100ng"];
cond3 = ["WT","D","0ng","1o10ng","1ng","5ng"];
cond4 = ["WT","D","0ng","1o10ng1","1o10ng2","1ng1","1ng2","5ng1","5ng2"];

% Experiment 1
X = qPCR_eff.exp1;
L = letters1;

for c = 1:size(cond1,2)
    qPCR.exp1.(['Lpp' char(cond1(c))]).data = [X.([L(c) '1']), X.([L(c) '2']), X.([L(c) '3'])];
    qPCR.exp1.(['Lpp' char(cond1(c))]).eff = [X.([L(c) '1' 'Eff']), X.([L(c) '2' 'Eff']), X.([L(c) '3' 'Eff'])];
    qPCR.exp1.(['Lpp' char(cond1(c))]).base = [X.([L(c) '1' 'base']), X.([L(c) '2' 'base']), X.([L(c) '3' 'base'])];

    if (c == 1)
        qPCR.exp1.(['idnT' char(cond1(c))]).data = [X.([L(c) '8']), X.([L(c) '10'])];
        qPCR.exp1.(['idnT' char(cond1(c))]).eff = [X.([L(c) '8' 'Eff']), X.([L(c) '10' 'Eff'])];
        qPCR.exp1.(['idnT' char(cond1(c))]).base = [X.([L(c) '8' 'base']), X.([L(c) '10' 'base'])];
    else
        qPCR.exp1.(['idnT' char(cond1(c))]).data = [X.([L(c) '7']), X.([L(c) '8']), X.([L(c) '10'])];
        qPCR.exp1.(['idnT' char(cond1(c))]).eff = [X.([L(c) '7' 'Eff']), X.([L(c) '8' 'Eff']), X.([L(c) '10' 'Eff'])];
        qPCR.exp1.(['idnT' char(cond1(c))]).base = [X.([L(c) '7' 'base']), X.([L(c) '8' 'base']), X.([L(c) '10' 'base'])];
    end

end

% Experiment 2
X = qPCR_eff.exp2;
L = letters2;

for c = 1:size(cond2,2)
    qPCR.exp2.(['Lpp' char(cond2(c)) '1']).data = [X.([L(c) '1']), X.([L(c) '2']), X.([L(c) '3'])];
    qPCR.exp2.(['Lpp' char(cond2(c)) '1']).eff = [X.([L(c) '1' 'Eff']), X.([L(c) '2' 'Eff']), X.([L(c) '3' 'Eff'])];
    qPCR.exp2.(['Lpp' char(cond2(c)) '1']).base = [X.([L(c) '1' 'base']), X.([L(c) '2' 'base']), X.([L(c) '3' 'base'])];
    qPCR.exp2.(['Lpp' char(cond2(c)) '2']).data = [X.([L(c) '10']), X.([L(c) '11']), X.([L(c) '12'])];
    qPCR.exp2.(['Lpp' char(cond2(c)) '2']).eff = [X.([L(c) '10' 'Eff']), X.([L(c) '11' 'Eff']), X.([L(c) '12' 'Eff'])];
    qPCR.exp2.(['Lpp' char(cond2(c)) '2']).base = [X.([L(c) '10' 'base']), X.([L(c) '11' 'base']), X.([L(c) '12' 'base'])];

    qPCR.exp2.(['idnT' char(cond2(c)) '1']).data = [X.([L(c) '7']), X.([L(c) '8']), X.([L(c) '9'])];
    qPCR.exp2.(['idnT' char(cond2(c)) '1']).eff = [X.([L(c) '7' 'Eff']), X.([L(c) '8' 'Eff']), X.([L(c) '9' 'Eff'])];
    qPCR.exp2.(['idnT' char(cond2(c)) '1']).base = [X.([L(c) '7' 'base']), X.([L(c) '8' 'base']), X.([L(c) '9' 'base'])];
    qPCR.exp2.(['idnT' char(cond2(c)) '2']).data = [X.(['H' num2str(2*c-1)]), X.(['H' num2str(2*c)])];
    qPCR.exp2.(['idnT' char(cond2(c)) '2']).eff = [X.(['H' num2str(2*c-1) 'Eff']), X.(['H' num2str(2*c) 'Eff'])];
    qPCR.exp2.(['idnT' char(cond2(c)) '2']).base = [X.(['H' num2str(2*c-1) 'base']), X.(['H' num2str(2*c) 'base'])];

end

% Experiment 3
X = qPCR_eff.exp3;
L = letters3;

for c = 1:size(cond3,2)
    if (c == 3)
        qPCR.exp3.(['Lpp' char(cond3(c))]).data = [X.([L(c) '1']), X.([L(c) '2'])];
        qPCR.exp3.(['Lpp' char(cond3(c))]).eff = [X.([L(c) '1' 'Eff']), X.([L(c) '2' 'Eff'])];
        qPCR.exp3.(['Lpp' char(cond3(c))]).base = [X.([L(c) '1' 'base']), X.([L(c) '2' 'base'])];
    else
        qPCR.exp3.(['Lpp' char(cond3(c))]).data = [X.([L(c) '1']), X.([L(c) '2']), X.([L(c) '3'])];
        qPCR.exp3.(['Lpp' char(cond3(c))]).eff = [X.([L(c) '1' 'Eff']), X.([L(c) '2' 'Eff']), X.([L(c) '3' 'Eff'])];
        qPCR.exp3.(['Lpp' char(cond3(c))]).base = [X.([L(c) '1' 'base']), X.([L(c) '2' 'base']), X.([L(c) '3' 'base'])];
    end

    qPCR.exp3.(['idnT' char(cond3(c))]).data = [X.([L(c) '4']), X.([L(c) '5']), X.([L(c) '6'])];
    qPCR.exp3.(['idnT' char(cond3(c))]).eff = [X.([L(c) '4' 'Eff']), X.([L(c) '5' 'Eff']), X.([L(c) '6' 'Eff'])];
    qPCR.exp3.(['idnT' char(cond3(c))]).base = [X.([L(c) '4' 'base']), X.([L(c) '5' 'base']), X.([L(c) '6' 'base'])];

end

% Experiment 4
X = qPCR_eff.exp4;

for c = 1:size(cond4,2)
    if(c == 1)
        qPCR.exp4.(['Lpp' char(cond4(c))]).data = [X.(['B'  num2str(c)]), X.(['C'  num2str(c)])];
        qPCR.exp4.(['Lpp' char(cond4(c))]).eff  = [X.(['B' num2str(c) 'Eff']), X.(['C' num2str(c) 'Eff'])];
        qPCR.exp4.(['Lpp' char(cond4(c))]).base  = [X.(['B' num2str(c) 'base']), X.(['C' num2str(c) 'base'])];
    else
        qPCR.exp4.(['Lpp' char(cond4(c))]).data = [X.(['A'  num2str(c)]), X.(['B'  num2str(c)]), X.(['C'  num2str(c)])];
        qPCR.exp4.(['Lpp' char(cond4(c))]).eff  = [X.(['A' num2str(c) 'Eff']), X.(['B' num2str(c) 'Eff']), X.(['C' num2str(c) 'Eff'])];
        qPCR.exp4.(['Lpp' char(cond4(c))]).base  = [X.(['A' num2str(c) 'base']), X.(['B' num2str(c) 'base']), X.(['C' num2str(c) 'base'])];
    end

    qPCR.exp4.(['idnT' char(cond4(c))]).data = [X.(['D'  num2str(c)]), X.(['E'  num2str(c)]), X.(['F'  num2str(c)])];
    qPCR.exp4.(['idnT' char(cond4(c))]).eff = [X.(['D' num2str(c) 'Eff']), X.(['E' num2str(c) 'Eff']), X.(['F' num2str(c) 'Eff'])];
    qPCR.exp4.(['idnT' char(cond4(c))]).base = [X.(['D' num2str(c) 'base']), X.(['E' num2str(c) 'base']), X.(['F' num2str(c) 'base'])];

end
%% Determine the minimum variation in efficiency within the calculated values
% save fluorescence statistics at the minimum efficiency value.
% calculate ideal threshold given all the mean fluorescence values.
% NOTE** search the datasets within the linear range.

conds.cond1 = cond1;
conds.cond2 = ["WT1","D1","0ng1","10ng1","50ng1","100ng1","WT2","D2","0ng2","10ng2","50ng2","100ng2"];
conds.cond3 = cond3;
conds.cond4 = cond4;

for exp = 1:4
    for c = 1:size(conds.(['cond' num2str(exp)]),2)
        experiment = (['exp' num2str(exp)]);
        gene1 = (['Lpp' char(conds.(['cond' num2str(exp)])(c))]);
        gene2 = (['idnT' char(conds.(['cond' num2str(exp)])(c))]);

        qPCR.(experiment).(gene1).effWinds = window_calc(qPCR.(experiment).(gene1).data);
        qPCR.(experiment).(gene2).effWinds = window_calc(qPCR.(experiment).(gene2).data);

        [~,I1] = min(qPCR.(experiment).(gene1).effWinds(:,5));
        [~,I2] = min(qPCR.(experiment).(gene2).effWinds(:,5));

        qPCR.(experiment).(gene1).effStats = [qPCR.(experiment).(gene1).effWinds(I1,1), qPCR.(experiment).(gene1).effWinds(I1,3)];
        qPCR.(experiment).(gene2).effStats = [qPCR.(experiment).(gene2).effWinds(I2,1), qPCR.(experiment).(gene2).effWinds(I2,3)];

    end
end

%% Interpolate the corrected fluorescence data, and determine when fluorescence crosses the threshold(Ct)
% For each experiment, the threshold is the lowest overall value from
% plugging in each efficiency index.
% 

%% Calculating threshold
for exp = 1:4
    minF = 1000; % Just setting it to something obviously too high for the data in these experiments
    for c = 1:size(conds.(['cond' num2str(exp)]),2)
        %
        experiment = (['exp' num2str(exp)]);
        gene1 = (['Lpp' char(conds.(['cond' num2str(exp)])(c))]);
        gene2 = (['idnT' char(conds.(['cond' num2str(exp)])(c))]);

        F_Lpp = qPCR.(experiment).(gene1).effStats(1);
        F_idnT = qPCR.(experiment).(gene2).effStats(1);

        tempF = [F_Lpp, F_idnT];
        tempMin = min(tempF);

        if(tempMin < minF)
            minF = tempMin;
        end
    end
    F_thresh(exp) = minF;
end

% Interpolating Ct values

for exp = 1:4
    for c = 1:size(conds.(['cond' num2str(exp)]),2)
        %
        experiment = (['exp' num2str(exp)]);
        gene1 = (['Lpp' char(conds.(['cond' num2str(exp)])(c))]);
        gene2 = (['idnT' char(conds.(['cond' num2str(exp)])(c))]);

        %% This section needs help.
        I_lpp = 1;
        I_idnT = 1;

        for k = 1:40
            if(qPCR.(experiment).(gene1).data(k,1) < F_thresh(exp))
                I_lpp = I_lpp + 1;
            end

            if(qPCR.(experiment).(gene2).data(k,1) < F_thresh(exp))
                I_idnT = I_idnT + 1;
            end
        end
%%
        
        % Lpp
        for i = 1:size(qPCR.(experiment).(gene1).data,2)
            qPCR.(experiment).(gene1).cycle(i) = interp1(qPCR.(experiment).(gene1).data((I_lpp-2):(I_lpp+2),i),(I_lpp-2):(I_lpp+2),F_thresh(exp),'spline');
        end
        % idnT
        for i = 1:size(qPCR.(experiment).(gene2).data,2)
            qPCR.(experiment).(gene2).cycle(i) = interp1(qPCR.(experiment).(gene2).data((I_idnT-2):(I_idnT+2),i),(I_idnT-2):(I_idnT+2),F_thresh(exp),'spline');
        end
    end
end

%% Using the Ct values and efficiencies, calculate fold changes for each condition relative to WT
for exp = 1:4
    for c = 1:size(conds.(['cond' num2str(exp)]),2)
        %
        experiment = (['exp' num2str(exp)]);
        gene1 = (['Lpp' char(conds.(['cond' num2str(exp)])(c))]);
        gene2 = (['idnT' char(conds.(['cond' num2str(exp)])(c))]);

        if (exp == 2)
            if (c < 7)
                Eff_LppWT = qPCR.(experiment).LppWT1.effStats(2);
                C_LppWT = mean(qPCR.(experiment).LppWT1.cycle);
                Eff_idnTWT = qPCR.(experiment).idnTWT1.effStats(2);
                C_idnTWT = mean(qPCR.(experiment).idnTWT1.cycle);
            else
                Eff_LppWT = qPCR.(experiment).LppWT2.effStats(2);
                C_LppWT = mean(qPCR.(experiment).LppWT2.cycle);
                Eff_idnTWT = qPCR.(experiment).idnTWT2.effStats(2);
                C_idnTWT = mean(qPCR.(experiment).idnTWT2.cycle);
            end
        else
            Eff_LppWT = qPCR.(experiment).LppWT.effStats(2);
            C_LppWT = mean(qPCR.(experiment).LppWT.cycle);
            Eff_idnTWT = qPCR.(experiment).idnTWT.effStats(2);
            C_idnTWT = mean(qPCR.(experiment).idnTWT.cycle);
        end

        Eff_LppCond = qPCR.(experiment).(gene1).effStats(2);
        C_LppCond = mean(qPCR.(experiment).(gene1).cycle);
        Eff_idnTCond = qPCR.(experiment).(gene2).effStats(2);
        C_idnTCond = mean(qPCR.(experiment).(gene2).cycle);

        RatioWT = (Eff_LppWT^C_LppWT)/(Eff_idnTWT^C_idnTWT);
        RatioCond = (Eff_idnTCond^C_idnTCond)/(Eff_LppCond^C_LppCond);

        qPCR.(experiment).(gene1).FC = RatioCond*RatioWT;
    end
end



%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%



%% Functions
function data = qPCR_load(FILENAME,sheet,letters,list)
% fct start

% These define the range of cells where the data is loaded from.
lower = 4;
upper = 43;

for i = 1:size(letters,2)
    for j = 1:size(list,2)
        suffix = [letters(i) num2str(list(j))];
        cells = ['C' num2str(lower) ':' 'C' num2str(upper)];

        data.(suffix) = cell2mat(readcell(FILENAME,"Sheet",sheet,"Range",cells));

        lower = lower + 42;
        upper = upper + 42;
    end
end

% fct end
end
%%
function [data, baseline] = BL_Correct(X)
% fct start

X = X/200; %% Data normalized to be in the range the algorithm in ref [1] was developed for. Units are arbitrary.

% Find ranges for calculating slopes
[imax, imin, i_0] = deriv_2(X); 

% Guess baseline
guess = min(X); 
dif_slope = BL_iterator(X,guess,imin,i_0,imax);

while(dif_slope < 0)
    guess = guess*0.99;
    dif_slope = BL_iterator(X,guess,imin,i_0,imax);
end

step = 0.0005*guess;

eps = abs(dif_slope); 

%
while(eps > 0.00005)
    if(step > 0)
        if(dif_slope > 0)
            guess = guess + step;
            dif_slope = BL_iterator(X,guess,imin,i_0,imax);
            eps = abs(dif_slope);
        elseif(dif_slope < 0)
            guess = guess - 2*step;
            step = step*0.5;
            dif_slope = BL_iterator(X,guess,imin,i_0,imax);
            eps = abs(dif_slope);
        else
            eps = 0;
            %guess = min(X);
        end
    else
        eps = 0;
    end
    %eps
    %guess
end
%
baseline = guess;

data = X - baseline;
% fct end
end

%%
function dif_slope = BL_iterator(X, guess,imin,i_0,imax)
%
XBC = X - guess;

log_XBC = real(log10(XBC));

%
i_nn = imax; %nn = not negative
i_nnTemp = i_nn;

for i = 1:8
    if(XBC((imax+1)-i) > 0)
        i_nnTemp = (imax+1-i);
    end
end


%%
iBot = imax-5;

Y_low = (log_XBC((iBot):(iBot+3)) - log_XBC(iBot));
Y_high = (log_XBC((iBot+3):(iBot+6)) - log_XBC(iBot+3));

X_low = transpose((iBot):(iBot + 3)) - (iBot);
X_high = transpose((iBot + 3):(iBot+6)) - (iBot + 3);

slope_low = X_low\Y_low;
slope_high = X_high\Y_high;


dif_slope = (slope_high - slope_low)/(slope_high);

%
end

%%
function [imax,imin,i_0] = deriv_2(X)
%
deriv = zeros(39,1);
% calculate second derivative
for i = 2:39
    deriv(i) = X(i-1) - 2*X(i) + X(i+1);
end
% Find the min and max of the second derivative to determine the linear
% region
[m1, imax] = max(deriv);
[m2, imin] = min(deriv);

i_0 = ceil((imax+imin)/2);
%
end
%%
function data = Baseline_calc(letters, list, qPCR)
%
for i = 1:size(letters,2)
    for j = 1:size(list,2)
        suffix = [letters(i) num2str(list(j))];

        [data.(suffix), data.([suffix 'base'])] = BL_Correct(qPCR.(suffix));
        %suffix;
    end
end
%
end
%%
function data = Efficiencies(letters,list,qPCR)
%
for i = 1:size(letters,2)
    for j = 1:size(list,2)
        suffix = [letters(i) num2str(list(j))];

        data.([suffix 'Eff']) = Efficiency_calc(qPCR.(suffix));
        data.(suffix) = qPCR.(suffix);
        data.([suffix 'base']) = qPCR.([suffix 'base']);
        %suffix;
    end
end
%
end
%%
function data = Efficiency_calc(qPCR)
%
data = zeros(40,1);
data(1) = 1;

for i = 2:40
    data(i) = qPCR(i)/qPCR(i-1);
end
%
end
%%
function data = window_calc(inputs)
%
Nreps = size(inputs,2);
imax = zeros(Nreps,1);
imin = zeros(Nreps,1);
i_0 = zeros(Nreps,1);
Eff = zeros(40,Nreps);
F_avg = 0;
%F_low = [];
Eff_avg = 0;

% Initial window calculation
for i = 1:Nreps
    [imax(i),imin(i),i_0(i)] = deriv_2(inputs(:,i));
    Eff(:,i) = Efficiency_calc(inputs(:,i));
    F_avg = F_avg + inputs(imax(i)+2,i)/Nreps;
    Eff_avg = Eff_avg + Eff(imax(i),i)/Nreps;
    %F_low = [F_low,inputs(imax(i)-4,i)];
end
upper = F_avg;
lower = upper - (Eff_avg)^4;
%lower = min(F_low);
% Call function to generate a new window
upper_old = 0;
data = [];
mean_eff = 1;
std_eff = 1;

while(upper ~= upper_old)
    data = [data;upper,lower,mean_eff,std_eff,(std_eff/mean_eff)];
    [upper,lower,upper_old,mean_eff,std_eff] = new_window(inputs,upper,lower,Eff,Nreps,imax);
    upper_old;
end
%
end
%%
function [upper,lower,upper_old,mean_eff,std_eff] = new_window(inputs,upper,lower,Eff,Nreps,imax)
%

Avg_eff = zeros(Nreps,1);
for R = 1:Nreps
    counter = 0;
    temp_eff = 0;
    for i = (imax(R)-4):40
        if(inputs(i,R)>=lower)
            if(inputs(i,R)<=upper)
                counter = counter+1;
                temp_eff = temp_eff + Eff(i,R);
            end
        end
    end
    counter;
    if(counter>0)
        Avg_eff(R) = temp_eff/counter;
        exit = false;
    else
        exit = true;
        R = Nreps;
    end
    exit;
end
mean_eff = mean(Avg_eff);
std_eff = std(Avg_eff);

if(exit == false)
    upper_old = upper;
%
    upper = upper - mean_eff/2;
    lower = upper - (mean_eff)^4;

else
    upper_old = upper;
end

[upper,lower,upper_old,mean_eff,std_eff];
%
end
%%
function [cycle, mean_eff, CV_eff] = Eff_stats(qPCR)
%[cycle, mean_eff, CV_eff, twoCycles]

[imax,imin,i_0] = deriv_2(qPCR.data(:,2));

mean_eff = mean(mean(qPCR.eff( (imax-4:imax) , :)));
CV_eff = std(mean(qPCR.eff( (imax-4:imax) , :)));
cycle = imax;


end


