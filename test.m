N = 500 ;%total population
I0 = 100; % initial number of infected
T = 10; % period 
beta=zeros(1,T);
gamma=zeros(1,T);

beta(1)= (1.1)*10^-9; %infection rate
gamma(1) = 0.0101; %recovery rate

dt = 1; % time interval of 6 hours 

% Calculate the model
[S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak

tt = 0:dt:T-dt;

% Curve
plot(tt,S,'b--o',tt,I,'r--o',tt,R,'g--o','LineWidth',2); 
grid on;
xlabel('Days');
ylabel('Number of individuals');

legend('S','I','R');

for jj=1:T
fprintf('Value of parameter R0 of day %d is %.2f',jj,N.*beta(jj)./gamma(jj))
disp(' ')
end


function [S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt)


%random numbers

    I = [50,75,100,125,109,100,99,70,50,40];
     S = zeros(1,T/dt);
    S(1) = N-I(1);

    R = [0,12,13,30,49,55,80,110,140,165];

     for tt = 1:T-1

        min=10;

       ans=fsolve(@(x) R(tt+1)-(x.*I(tt)).*dt-R(tt),[0 0.5 1 1.5 2 2.5 3]);
      for ii=1:7
        if ans(ii)>0 &&ans(ii)<4         
               gamma(tt+1)=ans(ii);
        end
       end

       min=10;

       ans=fsolve(@(y) I(tt+1)-(y.*I(tt).*(N-I(tt)-R(tt)) - gamma(tt+1).*I(tt)).* dt-I(tt),[0 0.5 1 1.5 2 2.5 3]);
       for ii=1:7
        if ans(ii)>0 &&ans(ii)<4
               beta(tt+1)=ans(ii);
           end
       end
        R(tt+1)=R(tt)+(gamma(tt+1).*I(tt)).*dt
        I(tt+1)=(beta(tt+1).*I(tt).*(N-I(tt)-R(tt)) - gamma(tt+1).*I(tt)).* dt+I(tt)
        S(tt+1) = N -I(tt+1)-R(tt+1)
        %check values

     end
end
