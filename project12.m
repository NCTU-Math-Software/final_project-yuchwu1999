
N = 23563356 ;%total population
I0 = 0; % initial number of infected
T = 331; % period
beta=zeros(1,T);
gamma=zeros(1,T);

beta(1)= 1.*10^(-9); %infection rate
gamma(1) =0.01; %recovery rate

dt = 1; % time interval of 6 hours

% Calculate the model
[S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak

tt = 0:dt:T-dt;

% Curve
plot(tt,I,'r',tt,R,'g','LineWidth',2);

grid on;
xlabel('Days');
ylabel('Number of individuals');

legend('I','R');

for jj=1:T
fprintf('Value of parameter R0 of day %d is %.2f',jj,(N.*beta(jj))./gamma(jj))
disp(' ')
end


function [S,I,R,beta,gamma] = sir_model(beta,gamma,N,I0,T,dt)




     S = zeros(1,T/dt);
%number of recovered people
    R = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5, 5, 5, 6, 9, 9, 12, 12, 12, 12, 12, 12, 13, 13, 15, 15, 17, 17, 20, 20, 20, 20, 22, 22, 22, 26, 28, 28, 28, 29, 29, 29, 29, 30, 30, 39, 39, 39, 45, 50, 50, 50, 54, 57, 61, 67, 80, 91, 99, 109, 114, 124, 137, 155, 166, 178, 189, 203, 217, 236, 253, 264, 275, 281, 290, 307, 311, 322, 324, 324, 332, 334, 334, 339, 347, 355, 361, 366, 368, 372, 375, 383, 387, 389, 395, 398, 401, 402, 407, 408, 411, 414, 415, 416, 419, 420, 420, 421, 423, 427, 427, 428, 428, 429, 429, 430, 430, 430, 431, 431, 431, 431, 431, 433, 433, 434, 434, 434, 434, 435, 435, 435, 435, 435, 435, 435, 435, 435, 437, 437, 438, 438, 438, 438, 438, 438, 438, 438, 438, 438, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 441, 441, 441, 441, 441, 441, 441, 443, 443, 443, 443, 443, 443, 443, 450, 450, 450, 450, 450, 450, 450, 457, 457, 457, 457, 457, 457, 457, 462, 462, 462, 462, 462, 462, 462, 471, 471, 471, 473, 473, 475, 475, 475, 475, 475, 475, 475, 476, 476, 477, 478, 478, 479, 479, 479, 479, 480, 480, 480, 480, 480, 482, 482, 483, 484, 484, 484, 484, 485, 485, 486, 487, 488, 488, 488, 489, 489, 491, 491, 491, 491, 491, 493, 495, 495, 497, 497, 502, 502, 502, 502, 508, 513, 514, 515, 518, 519, 521, 521, 523, 523, 523, 524, 526, 528, 528, 532, 533, 535, 536, 536, 539, 541, 545, 546, 546, 548, 549, 549, 553, 555, 555, 556, 565, 565,565, 568, 570, 572, 572, 574, 574, 574, 582, 585, 590, 595, 601, 606, 606];

    %number of infected people
    I=R-[0, 1, 1, 1, 3, 3, 4, 5, 8, 8, 9, 10, 10, 10, 10, 11, 11, 16, 16, 17, 18, 18, 18, 18, 18, 18, 18, 20, 22, 22, 24, 24, 26, 26, 28, 30, 31, 32, 32, 34, 39, 40, 41, 42, 42, 44, 45, 45, 45, 45, 47, 48, 49, 50, 53, 59, 67, 77, 100, 108, 135, 153, 169, 195, 216, 235, 252, 267, 283, 298, 306, 322, 329, 339, 348, 355, 363, 373, 376, 379, 380, 382, 385, 388, 393, 393, 395, 395, 395, 398, 420, 422, 425, 426, 427, 428, 429, 429, 429, 429, 429, 429, 429, 432, 436, 438, 438, 439, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 441, 441, 441, 441, 441, 441, 441, 441, 442, 442, 442, 443, 443, 443, 443, 443, 443, 443, 443, 443, 443, 443, 443, 443, 443, 445, 445, 445, 446, 446, 446, 446, 446, 446, 446, 447, 447, 447, 447, 447, 447, 447, 448, 449, 449, 449, 449, 449, 449, 449, 451, 451, 451, 451, 451, 451, 452, 454, 454, 455, 455, 455, 455, 455, 458, 458, 458, 462, 467, 467, 467, 467, 474, 475, 475, 476, 476, 477, 477, 479, 480, 480, 480, 481, 481, 481, 482, 484, 485, 486, 486, 486, 487, 487, 487, 487, 487, 487, 487, 487, 488, 488, 488, 488, 489, 489, 490, 492, 493, 494, 495, 495, 496, 498, 498, 498, 499, 499, 500, 503, 503, 506, 507, 509, 509, 509, 509, 510, 510, 510, 513, 513, 514, 515, 517, 517, 517, 518, 521, 523, 524, 527, 527, 527, 529, 529, 529, 530, 534, 534, 534, 539, 542, 543, 547, 547, 549, 549, 549, 549, 550, 553, 554, 555, 558, 563, 567, 568, 569, 573, 573, 577, 578, 580, 584, 589, 597, 600, 602, 603, 605, 607, 609, 611, 611, 617, 618, 618, 623, 625, 639, 648, 651, 675,675, 679, 685, 686, 690, 694, 716, 716, 718, 720, 724, 725, 733, 736, 740];
    I=I.*(-1);
    S(1) = N-I(1);
    for tt = 1:T-1

       gamma(tt+1)=fsolve(@(x) (R(tt+1)-(x*I(tt))-R(tt) * dt),1);
       beta(tt+1)=fsolve(@(y) (I(tt+1)-(y*I(tt)*(N-I(tt)-R(tt)) - gamma(tt)*I(tt)) * dt-I(tt)),1);



        S(tt+1) = N -I(tt)-R(tt);




    end
end
