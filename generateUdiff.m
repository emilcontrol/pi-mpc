% This script plots the surface plot with e and sum(e) on x-, and y-axis
% and difference between u_pi and u_mpc on the z-axis.

clear all
[ad,bd,A,B_u,B_d,R,Ts,N,Ndist,K,Ti,Kp,Ki,K_fb,distVec,distVec_ctrl,x0] = readParamsToWorkspace();

% Defining the grid of e, sum(e)
n = 101;
xminE = -10;
xmaxE = 10;
xminIE = -1000;
xmaxIE = 1000;

errGrid = linspace(xminE,xmaxE,n);
IntErrGrid = linspace(xminIE,xmaxIE,n);

% Data vectors for plot
uPI = zeros(n,n,3);
uMPC = zeros(n,n,3);
uDiff = zeros(n,n,3);

% data to overleaf. col1:e, col2:sum(e), col3:u_pi-u_mpc
dataToOverleaf = zeros(n*n,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stuff for MPC
[Q, Qf] = calculateQdiscForPI(ad,bd,Ki,Kp,R);

controlHorizon = 60;
predHorizon = controlHorizon;

% Optimization matrices H and fnx. G and F we need for constraint matrices
[H,fnx_u,fnx_d, G_u, G_d,F] = CalculateQPMtx(A,B_u,B_d,Q,R,Qf,controlHorizon,predHorizon);

% Hard constraints only for MPC
Xbounds = [-inf  inf ; % x1min, x1max
           -inf  inf]; % x2min, x2max
Ubounds = [-inf  inf]; % umin, umax

% Constraint matrices.
[Cbar, cnx, Fext, Gdext] = CalculateConstraintMtx(A, B_u, ...
                            Xbounds, Ubounds, controlHorizon, predHorizon, G_u, G_d, F);

% Suppress outputs from quadprog
OptOptions = optimoptions('quadprog', 'Display', 'off','OptimalityTolerance',eps, ...
    'MaxIterations',1000);
warning('off', 'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


overleafIndex = 1;
D = zeros(predHorizon,1);

for i = 1:n
    for j = 1:n
        % PI
        uPI(i,j,1) = errGrid(i);
        uPI(i,j,2) = IntErrGrid(j);
        uPI(i,j,3) = Kp*errGrid(i) + Ki*IntErrGrid(j);

        U = SolveMPC([IntErrGrid(j);errGrid(i)],H,fnx_u,fnx_d, Cbar,cnx,Fext,Gdext,D,zeros(controlHorizon,1), OptOptions);

        uMPC(i,j,1) = errGrid(i);
        uMPC(i,j,2) = IntErrGrid(j);
        uMPC(i,j,3) = U(1,:);

        % Diff
        uDiff(i,j,1) = errGrid(i);
        uDiff(i,j,2) = IntErrGrid(j);
        uDiff(i,j,3) = uPI(i,j,3) - uMPC(i,j,3);

        % Data format for overleaf plots
        dataToOverleaf(overleafIndex,1) = errGrid(i);
        dataToOverleaf(overleafIndex,2) = IntErrGrid(j);
        dataToOverleaf(overleafIndex,3) = uPI(i,j,3) - uMPC(i,j,3);
        overleafIndex = overleafIndex + 1;
    end
    disp(['Iteration: ', num2str(i), ' of ', num2str(n)]); 
end


% Surface plot of u_pi-u_mpc
surf(uDiff(:,:,1),uDiff(:,:,2),uDiff(:,:,3))
xlabel('Process Error'); ylabel('Integral of Process Error'); zlabel('Diff of uMPC and uPID');

% Some kind of evaluation number
disp(['DiffMtx: ',num2str(sum(sum(uDiff)))])

%%
% Save data as txt file
T_cell = table(dataToOverleaf(:,1),dataToOverleaf(:,2),dataToOverleaf(:,3), 'VariableNames',["e","sum_e","u_diff"]);
writetable(T_cell,'txtData/Udiff.txt');