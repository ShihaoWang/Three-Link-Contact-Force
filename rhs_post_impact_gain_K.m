function zdot = rhs_post_impact_gain_K(t,z,p)

% x = [ Theta; Alpha; Beta; Gama];
% xdot = [ Thetadot; Alphadot; Betadot; Gamadot];

Theta = z(1);     Alpha = z(2);     Beta = z(3); 
Thetadot = z(4);  Alphadot = z(5);  Betadot = z(6);

K = p.K;   % Derisvative gain

A_fn = p.Post_A_fn;     %@(Alpha,Beta,Theta)
B_fn = p.Post_B_fn;     %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_beta)

g_beta_fn = p.g_beta_fn;   %@(Alpha,Beta)
f_beta_fn = p.f_beta_fn;   %@(Alpha,Beta,Betadot,Theta)

g_beta = g_beta_fn(Alpha,Beta);
f_beta = f_beta_fn(Alpha,Beta,Betadot,Theta);
u_beta = (-K * Betadot - f_beta)/g_beta;

A = A_fn(Alpha,Beta,Theta);

B = B_fn(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_beta);  

xdot = z(4:6,1);
stateddot = -A\B;
xddot = stateddot(1:3,1);
zdot = [xdot;xddot];

end

