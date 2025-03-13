function [H,fnx_u,fnx_d, G_u, G_d,F] = CalculateQPMtx(Phi,Gamma_u,Gamma_d,Q,R,Qf,controlHorizon,predHorizon)

xdim = size(Phi,1);
udim = size(Gamma_u,2);

% Extended cost matrices. Qf and Q in Qbar. R in Rbar
Qbar = zeros(predHorizon*xdim); % +xdim for Qf
Rbar = zeros(controlHorizon*udim);

% Simulation matrices.xhat = Fxhat + Gu
F = zeros(predHorizon*xdim,xdim);

% G ÄR FEL OM HP INTE ÄR HC!!! Kanske löst nu.
G_u = zeros(predHorizon*xdim,controlHorizon*udim);
G_d = zeros(predHorizon*xdim,controlHorizon*udim);

% Assign values into Qbar- and F-matrices
for k = 1:predHorizon
    Qbar((k-1)*xdim+1:k*xdim , (k-1)*xdim+1:k*xdim) = Q;

    F((k-1)*xdim+1:k*xdim,1:xdim) = Phi^k;
end
Qbar((predHorizon-1)*xdim+1:(predHorizon)*xdim , (predHorizon-1)*xdim+1:(predHorizon)*xdim) = Qf;

% Assign values into Rbar and G-matrices
for k = 1:controlHorizon
    Rbar((k-1)*udim+1:k*udim , (k-1)*udim+1:k*udim) = R;
    for j = 1:k
        G_u((k-1)*xdim+1:k*xdim , (j-1)*udim+1:j*udim) = Phi^(k-j)*Gamma_u;
        G_d((k-1)*xdim+1:k*xdim , (j-1)*udim+1:j*udim) = Phi^(k-j)*Gamma_d;
    end
end

% If predHorizon > controlHorizon, add more Phi^k*Gamma
for k = controlHorizon+1:predHorizon
    for j = 1:controlHorizon
        G_u((k-1)*xdim+1:k*xdim , (j-1)*udim+1:j*udim) = Phi^(k-j)*Gamma_u;
        G_d((k-1)*xdim+1:k*xdim , (j-1)*udim+1:j*udim) = Phi^(k-j)*Gamma_d;
    end
end

H = G_u'*Qbar*G_u + Rbar;
fnx_u = F'*Qbar*G_u;
fnx_d = G_d'*Qbar*G_u;






