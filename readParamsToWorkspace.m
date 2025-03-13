function [ad,bd,A,B_u,B_d,R,Ts,N,Ndist,K,Ti,Kp,Ki,K_fb,distVec,distVec_ctrl,x0] = readParamsToWorkspace()
% This function reads all relevant parameters for the simulations.
% If all of this is located in each script, it becomes messy I think.

% Continous time
a=-0.0218;
b=-0.0521;

% ZOH discretized first-order system.
% Note that xdot=-ax+bu in cont. time (note the minus sign)
Ts = 1;
ad = exp(a*Ts);
bd = (b/a)*(exp(a*Ts)-1);
v=3.54*10^(-6);

% Extended ss-matrices in discrete time [sum(e), e]'
A = [1 1;0 ad];
B_u = [0;-1*bd];
B_d = [0;-v];

% R is control energy weight
R=1;

% Initial state
x0=[0;5];

% Simulation steps
N=2000;
%N=1200;

% When a disturbance enters
Ndist=round(N/2);
%Ndist=200;

% Disturbance
distVec = zeros(2*N,1);
d0 = 60000;

% Predicted disturbances. Change here to get known disturbances!
% The known disturbances. 0 gives unknown dist and d0 known.
% Change d0_ctrl to 0 for non-measured unknown disturbances.
distVec_ctrl = zeros(2*N,1);
d0_ctrl = 0;% d0; 

% Fill disturbance vectors
for i = 1:2*N
    if i > Ndist
        distVec(i) = d0;
        distVec_ctrl(i) = d0_ctrl;
    end
end

% Control parameters with different parametrizations
K=-0.9; % Negative since high u lets water out.
Ti=85; 
Kp=K;
Ki=K*Ts/Ti;

K_fb = [-1*Ki, -1*Kp]; %K positiv. u=-K ger negativa K då ki och kp negativa
end