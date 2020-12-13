

N = (2.3).*10^7 ;%total population
I0 = 100; % initial number of infected
T=2000;%period
beta=zeros(1,T);
gamma=zeros(1,T);

beta(1)= (1.1)*10^-9; %infection rate
gamma(1) = 0.0101; %recovery rate

dt = 1; % time interval of 6 hours (1/8 of a day)

% Calculate the model
[S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak

tt = 0:dt:T-dt;

% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2); 
grid on;
xlabel('Days');
ylabel('Number of individuals');

legend('S','I','R');

for jj=1:T
fprintf('Value of parameter R0 is %.2f',(N.*beta(jj))./gamma(jj))
disp(' ')
end


function [S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt)
   
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    
    for tt = 1:T-1
        
        dS = (-beta(tt)*I(tt)*S(tt)) * dt;
        dI = (beta(tt)*I(tt)*S(tt) - gamma(tt)*I(tt)) * dt;
        dR = (gamma(tt)*I(tt)) * dt;
       
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
        
        beta(tt+1)=beta(tt);
        gamma(tt+1)=gamma(tt);
        
       
       
    end
end
