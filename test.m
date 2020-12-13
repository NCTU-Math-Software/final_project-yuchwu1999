beta = 5*10^-9; %infection rate
gamma = 0.12; %recovery rate
N = 6.*10^7 %total population
I0 = 100; % initial number of infected
T = 365; % period 

dt = 1/8; % time interval of 6 hours (1/8 of a day)

fprintf('Value of parameter R0 is %.2f',(N*beta)./gamma)

% Calculate the model
[S,I,R] = sir_model(beta,gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak

tt = 0:dt:T-dt;

% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2); 
grid on;
xlabel('Days');
ylabel('Number of individuals');
legend('S','I','R');




function [S,I,R] = sir_model(beta,gamma,N,I0,T,dt)
   
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    for tt = 1:(T/dt)-1
        
        dS = (-beta*I(tt)*S(tt)) * dt;
        dI = (beta*I(tt)*S(tt) - gamma*I(tt)) * dt;
        dR = (gamma*I(tt)) * dt;
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
    end
end
