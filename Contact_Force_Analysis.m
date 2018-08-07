function Contact_Force_Analysis()

% This function is used to characterize the contact force during the whole
% stabilization process

global ratio
ratio = 0.5;

load('alpha_fn.mat');       load('theta_fn.mat');          load('alphadot_fn.mat');        load('thetadot_fn.mat');
% Initial Condition for contact force analysis
Wall_x = 0.1;       Wall_y = 0.1;          
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

[t,z] = ode23(@rhs_post_impact_AB_gain_K,tspan,Init_Condition,options,p);

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
    B_fn = p.Post_B_fn;     %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)
    
    h_fn = p.h_fn;   %@(Alpha,Beta,Betadot,Theta)
    g_alpha_fn = p.g_alpha_fn;   %@(Alpha,Beta)
    g_beta_fn = p.g_beta_fn;   %@(Alpha,Beta)
    
    g_alpha = g_alpha_fn(Alpha,Beta);
    g_beta = g_beta_fn(Alpha,Beta);
    h = h_fn(Alpha,Beta,Betadot,Theta);
    
    u_rhs = -p.K * Betadot - h;
    u_alpha = (1 - ratio) * u_rhs/(g_alpha);
    u_beta = ratio * u_rhs/(g_beta);
    
    
    A = A_fn(Alpha,Beta,Theta);
    B = B_fn(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha, u_beta);
    stateddot = -A\B;
    
    Thetaddot = stateddot(1);
    Alphaddot = stateddot(2);
    Betaddot = stateddot(3);
    
    
    FAx_Part_i = FAx_Part(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_alpha,u_beta);
    FAy_Part_i = FAy_Part(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_alpha,u_beta);
    
    FDx_Part_i = FDx_Part(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_alpha,u_beta);
    FDy_Part_i = FDy_Part(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_alpha,u_beta);
    
    
    FAx_Part_Tot = [FAx_Part_Tot; FAx_Part_i];
    FAy_Part_Tot = [FAy_Part_Tot; FAy_Part_i];
    FDx_Part_Tot = [FDx_Part_Tot; FDx_Part_i];
    FDy_Part_Tot = [FDy_Part_Tot; FDy_Part_i];
    
end


% Contact Force Analysis
% (Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta);
% (Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta);

%% Contact Force Analysis
% figure
% plot(t, FAx_Part_Tot,'LineWidth',1.5);
% hold on 
% plot(t, FAy_Part_Tot,'LineWidth',1.5);
% hold on 
% plot(t, FDx_Part_Tot,'LineWidth',1.5);
% hold on 
% plot(t, FDy_Part_Tot,'LineWidth',1.5);
% legend('FAx', 'FAy', 'FDx', 'FDy')
% 
% figure
% plot(t, FAx_Part_Tot./FAy_Part_Tot,'LineWidth',1.5);
% hold
% plot(t, FDy_Part_Tot./FDx_Part_Tot,'LineWidth',1.5);
% legend('Point A', 'Point D')s


%% In this case, the internal force points from the foot contact point to the hand contact point (Wall_x, Wall_y)
% As a result, the normalized internal force can be written as
F_Internal = [Wall_x, Wall_y];
F_Internal = F_Internal/norm(F_Internal);

M = F_Internal(1);                 N = F_Internal(2);


for i = 1:m
    
    FAx_Part_Tot = [FAx_Part_Tot; FAx_Part_i];
    FAy_Part_Tot = [FAy_Part_Tot; FAy_Part_i];
    FDx_Part_Tot = [FDx_Part_Tot; FDx_Part_i];
    FDy_Part_Tot = [FDy_Part_Tot; FDy_Part_i];
    
    A = FAx_Part_Tot(i);            B = FAy_Part_Tot(i);
    C = FDx_Part_Tot(i);            D = FDy_Part_Tot(i);
    
    x_low = ceil(max(-B/M, D/N));

    
    
    
    
    
end




end

