%% Lpp data fitting script
clear

%

load('Lpp_XY.mat')

%%
% Initial parameter guesses
R0 = 8.48;
R1 = 12;
beta = 1;
delta = 1;
gamma = -1;

% Setting what parameters are used in the fit
% 1: only beta
% 2: only delta
% 3: beta, delta, and epsilon

flag = 2;

% Bootstrapping parameters
list = zeros(11,1);
iterations = 1;


% Defining data to be used for fitting
data = [];
data(:,1) = Y;
data(:,2) = Yerr;

% Initial guesses, upper and lower bounds for fitting ub(i) = lb(i) fixes that value
X0 = [R0, R1, beta, delta, gamma];
lb = [8.48, 12, 0.00001, 0.00001, -10];
ub = [8.48, 12, 1000, 1000, -0.001];

% Set any additional options if required
options = optimset('Display','on','MaxIter',5000,'MaxFunEvals',2500,'DiffMinChange',0.1);

tic
% Generating indices used for bootstrapping
for i = 1:iterations
    for j=1:11
        list(j) = j;%randi(11,1); %j;
    end

    % Function handle
    fun = @(X) residual_Fcn(X,data,list,X3,GCR_interp2,flag);

    % Perform the fitting
    [xOpt(i,:),resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(fun, X0, lb, ub, options);

    i
end
toc

% Fit parameters output to screen
%xOpt(1) 
%xOpt(3)


%% Residuals Function
function residual = residual_Fcn(X,data,list,X3,counts,flag)

R_temp = RC_KP_Lpp(X(1),X(2),X(3),X(4),X(4),flag);

R_pred = interp1(0:0.0001:1,R_temp,X3,"pchip");

%index = zeros(11,1);
FC_pred = zeros(11,1);
count_interp = zeros(11,1);

baseline = round(X(2)*100)+1;

for j = 1:11
    %index(j) = round(R_pred(j)*100) + 1;

    count_interp(j) = interp1(0:0.01:40,counts(:,j),R_pred(j));

    FC_pred(j) = count_interp(j)/counts(baseline,11); %counts(index(j),j)/counts(baseline,11);
end

FC_data = data(:,1);

% Standard deviation matrix
SD = data(:,2);

Res = zeros(size(list,1),1);
Weight = zeros(size(list,1),1);

for i = 1:size(list,1)
    a = list(i);

    Res(i) = FC_pred(a) - FC_data(a); % unsquared difference, lsqnonlin does the squaring
    Weight(i) = 1/SD(a);
end

% Apply weights to the residuals
residual = Res.*Weight;

end