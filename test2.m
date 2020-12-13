N = (2.3).*10^7 ;%total population
I0 = 100; % initial number of infected
T = 10; % period
beta=zeros(1,T);
gamma=zeros(1,T);

beta(1)= (1.1)*10^-9; %infection rate
gamma(1) = 0.0101; %recovery rate

dt = 1; % time interval 
% Calculate the model
[S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak

tt = 0:dt:T-dt;



for jj=1:T
fprintf('Value of parameter R0 of day %d is %.2f',jj,(N.*beta(jj))./gamma(jj))
disp(' ')
end


function [S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt)

    S = zeros(1,T/dt);
    S(1) = N;
    I = [100,105,106,108,109,110,111,114,107,105];


    R = [0,12,13,14,15,16,17,18,19,20];


    for tt = 1:T-1

        gamma(tt+1)=fsolve(@(x) R(tt+1)-(x*I(tt))-R(tt) * dt,1);
       beta(tt+1)=fsolve(@(y) I(tt+1)-(y*I(tt)*(N-I(tt)-R(tt)) - gamma(tt)*I(tt)) * dt-I(tt),1);

        dI = (beta(tt+1)*I(tt)*(N-I(tt)-R(tt)) - gamma(tt+1)*I(tt)) * dt;
        dR = (gamma(tt+1)*I(tt)) * dt;



        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
        S(tt+1) = N -I(tt)-R(tt);




    end
end
