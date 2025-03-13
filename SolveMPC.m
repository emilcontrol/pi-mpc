function U = SolveMPC(xt,H,fnx_u,fnx_d,Cbar,cnx,Fext,Gdext,D,U_mpc_old, OptOptions)
% Solve quadratic programming problem
U = quadprog(H,fnx_u'*xt+fnx_d'*D, Cbar,cnx-Fext*xt-Gdext*D, [],[],[],[],U_mpc_old, OptOptions);
%U = quadprog(H,fnx'*xt, [],[], [],[],[],[],[], OptOptions);
end