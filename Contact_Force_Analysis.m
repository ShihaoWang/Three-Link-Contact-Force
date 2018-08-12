function Contact_Force_Analysis()

% This function is used to characterize the contact force during the whole
% stabilization process

global ratio
ratio = 0.5;

load('alpha_fn.mat');       load('theta_fn.mat');          load('alphadot_fn.mat');        load('thetadot_fn.mat');
% Initial Condition for contact force analysis
Wall_x = 0.15;       Wall_y = 0.12;
Beta = pi/2;        Betadot = -2;

Alpha = alpha_fn(Beta,Wall_x,Wall_y);
Theta = theta_fn(Beta,Wall_x,Wall_y);
Alphadot = alphadot_fn(Alpha(1),Beta,Betadot);
Thetadot = thetadot_fn(Alpha(1),Beta,Betadot);

Init_Condition = [Theta(1); Alpha(1); Beta(1); Thetadot; Alphadot; Betadot];

p = Robot_Component_InertiaNLength();
p = Function_Attach(p);
p.K = 3.5;
p.Post_Impact_State = Init_Condition;

% Post-Impact integration
tspan = linspace(0,2,51);              % Assume that the system will be stabilize within 3s
options = odeset('AbsTol',1e-5,'RelTol',1e-5,'MassSingular','yes');

[t,z] = ode23(@rhs_post_impact_gain_K,tspan,Init_Condition,options,p);

[m,n] = size(z);
figure
Three_Link_Plot(z);

load('FAx_Part.mat');       load('FAy_Part.mat');       load('FDx_Part.mat');           load('FDy_Part.mat');

FAx_Part_Tot = [];          FAy_Part_Tot = [];          FDx_Part_Tot = [];          FDy_Part_Tot = [];
for i = 1:m
    
    state_i = z(i,:);
    Theta = state_i(1);             Alpha = state_i(2);             Beta = state_i(3);
    Thetadot = state_i(4);          Alphadot = state_i(5);          Betadot = state_i(6);
    
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
    stateddot = -A\B;
    
    Thetaddot = stateddot(1);
    Alphaddot = stateddot(2);
    Betaddot = stateddot(3);
    
    
    FAx_Part_i = FAx_Part(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
    FAy_Part_i = FAy_Part(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
    
    FDx_Part_i = FDx_Part(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
    FDy_Part_i = FDy_Part(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
    
    
    FAx_Part_Tot = [FAx_Part_Tot; FAx_Part_i];
    FAy_Part_Tot = [FAy_Part_Tot; FAy_Part_i];
    FDx_Part_Tot = [FDx_Part_Tot; FDx_Part_i];
    FDy_Part_Tot = [FDy_Part_Tot; FDy_Part_i];
    
end


% Contact Force Analysis
% (Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
% (Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta);

%% Contact Force Analysis
figure
plot(t, FAx_Part_Tot,'LineWidth',1.5);
hold on
plot(t, FAy_Part_Tot,'LineWidth',1.5);
hold on
plot(t, FDx_Part_Tot,'LineWidth',1.5);
hold on
plot(t, FDy_Part_Tot,'LineWidth',1.5);
legend('FAx', 'FAy', 'FDx', 'FDy')

figure
plot(t, FAx_Part_Tot./FAy_Part_Tot,'LineWidth',1.5);
hold
plot(t, FDy_Part_Tot./FDx_Part_Tot,'LineWidth',1.5);
legend('Point A', 'Point D')


%% In this case, the internal force points from the foot contact point to the hand contact point (Wall_x, Wall_y)
% As a result, the normalized internal force can be written as
F_Internal = [Wall_x, Wall_y];
F_Internal = F_Internal/norm(F_Internal);

M = F_Internal(1);                 N = F_Internal(2);

load('r_Pos_fn.mat');              load('r_Neg_fn.mat');
load('a_fn.mat');
load('b_fn.mat');
load('c_fn.mat');
load('d_fn.mat');

p.r_Pos_fn = r_Pos_fn;             p.r_Neg_fn = r_Neg_fn;
p.a_fn = a_fn;      p.b_fn = b_fn;      p.c_fn = c_fn;      p.d_fn = d_fn;

mu_A = [];
mu_D = [];
for i = 1:m
    
    FAx_Part_i = FAx_Part_Tot(i);       FAy_Part_i = FAy_Part_Tot(i);
    FDx_Part_i = FDx_Part_Tot(i);       FDy_Part_i = FDy_Part_Tot(i);
    
    Par_i = [FAx_Part_i, FAy_Part_i, FDx_Part_i, FDy_Part_i];
    
    [A, B, C, D, E, F, G, H] = Contact_Force_Symbolization(Par_i, Wall_x, Wall_y, pi/2);
    
%     X = fzero(@(x)C^2*F^2*G*x*abs(A + B*x) - A^2*D*G^2*abs(E + F*x) - C^2*E^2*H*abs(A + B*x) + B^2*C*G^2*x*abs(E + F*x) + D^2*F^2*G*x^3*abs(A + B*x) - D^2*E^2*H*x^2*abs(A + B*x) + B^2*C*H^2*x^3*abs(E + F*x) - A^2*D*H^2*x^2*abs(E + F*x) + A*B*C*G^2*abs(E + F*x) + C^2*E*F*G*abs(A + B*x) - 2*C*D*E^2*H*x*abs(A + B*x) - C^2*E*F*H*x*abs(A + B*x) - A*B*D*G^2*x*abs(E + F*x) - 2*A^2*D*G*H*x*abs(E + F*x) + 2*C*D*F^2*G*x^2*abs(A + B*x) + A*B*C*H^2*x^2*abs(E + F*x) + D^2*E*F*G*x^2*abs(A + B*x) - A*B*D*H^2*x^3*abs(E + F*x) - D^2*E*F*H*x^3*abs(A + B*x) + 2*B^2*C*G*H*x^2*abs(E + F*x) + 2*C*D*E*F*G*x*abs(A + B*x) + 2*A*B*C*G*H*x*abs(E + F*x) - 2*C*D*E*F*H*x^2*abs(A + B*x) - 2*A*B*D*G*H*x^2*abs(E + F*x),0);
    
    [x_tot, Feasi_tot, mu_A_tot, mu_D_tot] = Inner_Force_All(A, B, C, D, E, F, G, H);
%     figure
%     plot(x_tot, mu_A_tot,'LineWidth',1.5)
%     hold on
%     plot(x_tot, mu_D_tot,'LineWidth',1.5)
%     hold on
    tic
    [ mu_A_i, mu_D_i, x_i] = Inner_Force_Analytic(A, B, C, D, E, F, G, H, p);
    toc
    
%     tic
%     [ mu_A_i, mu_D_i, x_i] = Inner_Force_Optimization(A, B, C, D, E, F, G, H, p);
%     toc
    mu_A = [mu_A; mu_A_i];
    mu_D = [mu_D; mu_D_i];
    
end
figure 
plot(t, mu_A, 'LineWidth',1.5);
hold on
plot(t, mu_D, 'LineWidth',1.5);
legend('mu_A', 'mu_D')


end


