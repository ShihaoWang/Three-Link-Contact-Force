function Contact_Force_Analysis()

% This function is used to characterize the contact force during the whole
% stabilization process

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
p.K = 3;

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
    B_fn = p.Post_B_fn;     %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)
    
    
    g1_fn = p.g1_fn;   %@(Alpha,Beta)
    f1_fn = p.f1_fn;   %@(Alpha,Beta,Betadot,Theta)
    
    g1_fn = g1_fn(Alpha,Beta);
    f1_fn = f1_fn(Alpha,Beta,Betadot,Theta);
    u_beta = (-K * Betadot - f1_fn)/g1_fn;
    
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

end

