function zdot = rhs_post_impact_AB_gain_K(t,z,p)
global ratio

% x = [ Theta; Alpha; Beta; Gama];
% xdot = [ Thetadot; Alphadot; Betadot; Gamadot];

Theta = z(1);     Alpha = z(2);     Beta = z(3); 
Thetadot = z(4);  Alphadot = z(5);  Betadot = z(6);

Kd = p.K;   % Derisvative gain
Kp = 3;

Post_Impact_State = p.Post_Impact_State;

A_fn = p.Post_A_fn;     %@(Alpha,Beta,Theta)
B_fn = p.Post_B_fn;     %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)


h_fn = p.h_fn;   %@(Alpha,Beta,Betadot,Theta)
g_alpha_fn = p.g_alpha_fn;   %@(Alpha,Beta)
g_beta_fn = p.g_beta_fn;   %@(Alpha,Beta)

h = h_fn(Alpha,Beta,Betadot,Theta);
g_alpha = g_alpha_fn(Alpha,Beta);
g_beta = g_beta_fn(Alpha,Beta);

u_rhs = -Kd * Betadot - h;
u_alpha = (1 - ratio) * u_rhs/(g_alpha);
u_beta = ratio * u_rhs/(g_beta);

A = A_fn(Alpha,Beta,Theta);

B = B_fn(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha, u_beta);  

xdot = z(4:6,1);
stateddot = -A\B;
xddot = stateddot(1:3,1);
zdot = [xdot;xddot];

end

