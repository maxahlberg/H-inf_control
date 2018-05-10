% This is a Matlab file for designing H_infinity controller (assignment 3
% of SD2231)
clear all; clc; close all;
s=tf('s');

% systme parameters
m=22000;   %kg
j=700e3;   %kgm^2
c=40e3;    %Ns/m
k=2*300e3; %N/m
L=6;       %m
%% Frequencies for the sine waves
f1 = 1;     %Hz
f = f1;
f2 = 8;     %Hz

%% State space model for skyhook contorl
Ask=[0 1 0 0
    -2*k/m 0 0 0
    0 0 0 1
    0 0 -2*k*L^2/j 0];
Bsk=[0 0 0 0
    k/m k/m -1/m -1/m
    0 0 0 0
    -L*k/j L*k/j L/j -L/j];
Csk=[0 1 0 0
    0 0 0 1];
Dsk=zeros(2,4);

%% H_inf using linmod syntax

% ----------Change those parameters to change the actual model, the train!
M = 22000;  %kg
J = 700000; %kgm^2
c1 = 40000; %Ns/m
c2 = c1;    %Ns/m
k1 = 600000;%N/m
k2 = k1;    %N/m
L1 = 6;     %m
L2 = L1;    %m
% ---- Tuning coefficients for skyhook ---
cz = 1.52*10^5; %Controlling the bounce
cx = 6.32*10^6;  %Controlling the pitch 
% ----------------------------------------

%state space: The same as skyhook
%% State space model Passive system MODELS OF ALL REAL SYSTEMS

A = [0                    1              0                0
    -(k1+k2)/M     -(c1+c2)/M  -(L2*k2-L1*k1)/M      -(c2*L2-c1*L1)/M
    0                    0              0                1
    (k1*L1-k2*L2)/J  (c1*L1-c2*L2)/J   -(k1*L1^2+k2*L2^2)/J  -(c1*L1^2+c2*L2^2)/J];


B = [0                    0              0                 0
    k1/M               c1/M           k2/M               c2/M
    0                    0              0                 0
    (-k1*L1)/J       (-c1*L1)/J     (k2*L2)/J          (c2*L2)/J];

C = [1, 0, 0, 0
    0, 0, 1, 0];
D = zeros(2,4);
%% State space model Skyhook

As = [0                   1              0                0
    -(k1+k2)/M            0      -(L2*k2-L1*k1)/M         0
    0                     0              0                1
    (k1*L1-k2*L2)/J       0    -(k1*L1^2+k2*L2^2)/J       0];


Bs = [0                    0              0                 0
    k1/M                k2/M           -1/M              -1/M
    0                     0              0                 0
    (-k1*L1)/J       (k1*L2)/J          L1/J              -L2/J];

Cs = [1, 0, 0, 0
     0, 1, 0, 0
     0, 0, 1, 0
     0, 0, 0, 1];
Ds = zeros(4,4);
%% Eigenvalue analysis
eigSKY = eig(Ask);
eigSKY1 = imag(eigSKY(1));
eigSKY2 = imag(eigSKY(3));
%%
      
%Weighting functions

%For penalizing actuator force
Wa1=(0.00175*s+1)/(0.00025*s+1);
Wa2=Wa1;

%For penalizing bounce and pitch motions
eps=1;
%The frequency that the Hinf controller will try to reduce the most
wnb= eigSKY1;            %Find the right equation or value for wnb
wnchi= eigSKY2;          %Find the right equation or value for wnchi
s1b=-eps+1i*sqrt(wnb^2-eps^2);
s2b=-eps-1i*sqrt(wnb^2-eps^2);
s1chi=-eps+1i*sqrt(wnchi^2-eps^2);
s2chi=-eps-1i*sqrt(wnchi^2-eps^2);
%-------------------------------------------------------------
kb = 2.899591300000000e+03;
kchi = 3.386473300000000e+04;
%-------------------------------------------------------------
%kb=input('Enter the gain for Wb = '); 
%kchi=input('Enter the gain for Wchi = ');
Wb=(kb*s1b*s2b)/((s-s1b)*(s-s2b));
Wchi=(kchi*s1chi*s2chi)/((s-s1chi)*(s-s2chi));
%%  ----- Plotting the weight functions ----
% [H1, Wout1] = freqresp(1/Wb);
% H1 = squeeze(H1);
% H1 = abs(H1);
% [H2, Wout2] = freqresp(1/Wchi);
% H2 = squeeze(H2);
% H2 = abs(H2);
% figure
% loglog(Wout1,H1,Wout2,H2)
%figure
%bode(1/Wb)
%loglog(W',1/Wb)
%hold on
%bode(1/Wchi)
%%
%Extracting the extended model
[A_Pe,B_Pe,C_Pe,D_Pe] = linmod('Extended_model');% state space parameters of the extended system: Pe
Pe=ss(A_Pe,B_Pe,C_Pe,D_Pe);
%PeTF = ss2tf(Pe);
%%
% [SV, W]= sigma(Pe);
% figure
% loglog(W',SV(1,:),W',SV(4,:),'LineWidth',1.5) %Largest and smallest singular values over W
% xlabel('Frequency [rad/s]'); ylabel('Magnitude [dB]')
% title('Singular value plot of the passive system')

%%
%Calculating the controller
ncont = 2;%Number of control inputs
nmeas = 2;%Number of measured outputs provided to the controller
Pe=minreal(Pe);%This syntax cancels pole-zero pairs in transfer
%functions. The output system has minimal order and the same response
%characteristics as the original model.
[K,Pec,gamma,info]=hinfsyn(Pe,nmeas,ncont,'method','lmi'); % for working with the error
[Ainf, Binf, Cinf, Dinf]=ssdata(K);

%Now use the controller K in your simulation

%% Simulating
sim('PB_skyhook_model.slx')

%% Plotting
% Excitation one Step on 0.03 m

%figure
subplot(3,2,1)
plot(tout,Z,tout,Z3,tout,Z6,'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Bounce [m]')
title('Bounce of the 3 models with disturbance')
%saveas(gcf,'bounce_step_rob_kMinus.eps','epsc')

subplot(3,2,2)
%figure
plot(tout,X*(180/pi),tout,X3*(180/pi),tout,X6*(180/pi),'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Pitch [Degrees]')
title('Pitch of the 3 models with disturbance')
%saveas(gcf,'pitch_step_rob_kMinus.eps','epsc')

% Excitaton f = 1 Hz

subplot(3,2,3)
%figure
plot(tout,Z1,tout,Z4,tout,Z7,'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Bounce [m]')
title('Bounce of the 3 models f = 1 Hz')
%saveas(gcf,'bounce_f1_all.eps','epsc')

subplot(3,2,4)
%figure
plot(tout,X1*(180/pi),tout,X4*(180/pi),tout,X7*(180/pi),'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Pitch [Degrees]')
title('Pitch of the 3 models f = 1 Hz')
%saveas(gcf,'pitch_f1_all.eps','epsc')


% Excitation f = 8 Hz
subplot(3,2,5)
%figure
plot(tout,Z2,tout,Z5,tout,Z8,'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Bounce [m]')
title('Bounce of the 3 models f = 8 Hz')
%saveas(gcf,'bounce_f8_all.eps','epsc')

subplot(3,2,6)
%figure
plot(tout,X2*(180/pi),tout,X5*(180/pi),tout,X8*(180/pi),'LineWidth',1.5)
legend('Skyhook','Passive','H-inf')
xlabel('Time [s]');ylabel('Pitch [Degrees]')
title('Pitch of the 3 models f = 8 Hz')
%saveas(gcf,'pitch_f8_all.eps','epsc')

%% Force plots skyhook
fmax = [10000, 10000];
fmin = [-10000, -10000];

x = [0, 5];

figure
subplot(3,1,1)
%plot(tout,Fa,tout,Fb,tout,Fa3,tout,Fb3,x,fmax,'--r',x,fmin,'--r')
plot(tout,Fa3,tout,Fb3,x,fmax,'--r',x,fmin,'--r')
title('Actuator forces Step'); %legend('Fa1 Sky','Fa2 Sky','Fa1 H-inf','Fa2 H-inf')
legend('Fa1', 'Fa2')
xlabel('Time'); ylabel('Newton')

subplot(3,1,2)
%plot(tout,Fa1,tout,Fb1,tout,Fa4,tout,Fb4,x,fmax,'--r',x,fmin,'--r')
plot(tout,Fa4,tout,Fb4,x,fmax,'--r',x,fmin,'--r')
title('Actuator forces f = 1'); legend('Fa1','Fa2')%,'Fa1 H-inf','Fa2 H-inf')
xlabel('Time'); ylabel('Newton')

subplot(3,1,3)
%plot(tout,Fa2,tout,Fb2,tout,Fa5,tout,Fb5,x,fmax,'--r',x,fmin,'--r')
plot(tout,Fa5,tout,Fb5,x,fmax,'--r',x,fmin,'--r')
title('Actuator forces f = 8'); legend('Fa1','Fa2')%,'Fa1 H-inf','Fa2 H-inf')
xlabel('Time'); ylabel('Newton')
%saveas(gcf,'actuator_force_rob_kMinus.eps','epsc')



