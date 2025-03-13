function [Cbar, cnx, Fext, Gdext] = CalculateConstraintMtx(Phi, Gamma_u, ...
                            Xbounds, Ubounds, controlHorizon, predHorizon, G_u, G_d, F)

% Constraints are written as Cbar*U <= c
% Cbar är [-G ; G ; -I ; I]
% c är [-x_constraints_min + F*x0 + G_d*D ; x_constraints_max - F*x0 -G_d*D; 
%       -u_constraints_min ; u_constraints_max]

xdim = size(Phi,1);
udim = size(Gamma_u,2);

% Create a representation of the state vector that makes sense:
% cvec måste vara [x1min(1), x2min(1), x1min(2), x2min(2), ... 
%              ... , x1min(predHorizon), x2min(predHorizon), ...
%               x1max(1), x2max(1), x1max(2), x2max(2), ... 
%              ... , x1max(predHorizon), x2max(predHorizon), ...
%              umin(1), umin(2), ... , umin(controlHorizon), ...
%              umax(1), umax(2), ... , umax(controlHorizon)]

cnx = zeros(predHorizon*xdim*2 + controlHorizon*udim*2,1);

for i = 1:predHorizon
    for j = 1:xdim
        % Min constraints for the xdim states
        cnx((i-1)*xdim+j, 1) = -Xbounds(j,1);

        % Max constraints for the xdim states
        cnx(predHorizon*xdim + (i-1)*xdim+j, 1) = Xbounds(j,2);
    end
end

for i = 1:controlHorizon
    for j = 1:udim
        % Min constraints for u
        cnx(predHorizon*xdim*2 + (i-1)*udim+j,1) = -Ubounds(j,1);

        % Max constraints for u
        cnx(predHorizon*xdim*2 + controlHorizon*udim + ...
            (i-1)*udim+j,1) = Ubounds(j,2);
    end
end

% Create the extended F-vector and fill it. (the mtx that shall be *x0)
% Shall be multiplied with x0 later.
Fext = zeros(predHorizon*xdim*2 + controlHorizon*udim*2, xdim);
Fext(1:predHorizon*xdim, 1:xdim) = -F;
Fext(predHorizon*xdim+1:2*predHorizon*xdim, 1:xdim) = F;

% Create the extended G_d vector for the constraints and fill it
Gdext = zeros(predHorizon*xdim*2 + controlHorizon*udim*2, predHorizon);
Gdext(1:predHorizon*xdim, 1:predHorizon) = -G_d;
Gdext(predHorizon*xdim+1:2*predHorizon*xdim, 1:predHorizon) = G_d;


% Create Cbar
Cbar = zeros(predHorizon*xdim*2 + controlHorizon*udim*2,controlHorizon*udim);

% Set G matrices into Cbar
Cbar(1:predHorizon*xdim, 1:controlHorizon*udim) = -G_u;
Cbar(predHorizon*xdim+1:2*predHorizon*xdim, 1:controlHorizon*udim) = G_u;

% Set Identity matrices into Cbar
Cbar(2*predHorizon*xdim+1:2*predHorizon*xdim+controlHorizon*udim, ...
     1:controlHorizon*udim) = -eye(controlHorizon*udim);

Cbar(2*predHorizon*xdim+controlHorizon*udim+1:2*predHorizon*xdim+2*controlHorizon*udim, ...
     1:controlHorizon*udim) = eye(controlHorizon*udim);
