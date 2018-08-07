function zdot = rhs_post_impact_gain_K(t,z,p)

% x = [ Theta; Alpha; Beta; Gama];
% xdot = [ Thetadot; Alphadot; Betadot; Gamadot];

Theta = z(1);     Alpha = z(2);     Beta = z(3); 
Thetadot = z(4);  Alphadot = z(5);  Betadot = z(6);

K = p.K;   % Derisvative gain

A_fn = p.Post_A_fn;     %@(Alpha,Beta,Theta)
B_fn = p.Post_B_fn;     %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)

    
g1_fn = p.g1_fn;   %@(Alpha,Beta)
f1_fn = p.f1_fn;   %@(Alpha,Beta,Betadot,Theta)

g1_fn = g1_fn(Alpha,Beta);
f1_fn = f1_fn(Alpha,Beta,Betadot,Theta);
u_beta = (-K * Betadot - f1_fn)/g1_fn;

A = A_fn(Alpha,Beta,Theta);

B = B_fn(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_beta);  

xdot = z(4:6,1);
stateddot = -A\B;
xddot = stateddot(1:3,1);
zdot = [xdot;xddot];

end

