function R = RC_KP_Lpp(RC0,RC1,beta,delta,gamma,flag)

X = 0:0.0001:1;

%Xinterp = [0,0.045696956,0.118508815,0.139855277,0.412953654,0.446415287,0.571090583,0.65421360,0.6647,0.561892571,1];


%% Rigidity scales like 1-exp(-beta*Lpp density)

% 1: only rigidity
% 2: only Pressure
% 3: Both rigidity and pressure

if flag == 1
    denom = 1-exp(-1*beta);
    for i = 1:size(X,2)
        numerator = 1 - exp(-beta*X(i));

        B = numerator/denom;

        R(i) = ((RC0^3)*(1-B)+(RC1^3)*(B))^(1/3);
    end

    %K = (R.^3)/(R(end)^3);

elseif flag == 2
    denom = 1-exp(-1*delta);
    for i = 1:size(X,2)
        numerator = 1 - exp(-delta*X(i));

        D = numerator/denom;

        R(i) = ((1-D)/(RC0^3) + D/(RC1^3))^(-1/3);
    end

    %P = (R(end)^3)/(R.^3);
elseif flag == 3

    B_denom = 1-exp(-1*beta);
    D_denom = 1-exp(-1*delta);
    for i = size(X,2)
        B_num = 1 - exp(-beta*X(i));
        D_num = 1 - exp(-delta*X(i));

        B = B_num/B_denom;
        D = D_num/D_denom;

        top = (RC0^3)*(1-B) + (RC1^3)*B*(1+gamma*D_num);
        bottom = (1 - D) + (1+ gamma*D_num)*D;

        R(i) = (top/bottom)^(1/3);

    end

end

end