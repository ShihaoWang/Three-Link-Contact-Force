function p = Function_Attach(p)

% This function will loads all the necessary functions to structure p

Pre_A = load('Pre_A_fn.mat'); %@(Alpha,Beta,Theta)
Pre_B = load('Pre_B_fn.mat'); %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)
p.Pre_A_fn = Pre_A.A_fn;                   
p.Pre_B_fn = Pre_B.B_fn;

Post_A = load('Post_A_fn.mat'); %@(Alpha,Beta,Theta)
Post_B = load('Post_B_fn.mat'); %@(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)
p.Post_A_fn = Post_A.A_fn;                 
p.Post_B_fn = Post_B.B_fn;

load('FAx_Part.mat')  %@(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta)
load('FAy_Part.mat')  %@(Alpha,Alphadot,Alphaddot,Beta,Betaddot,Theta,Thetadot,Thetaddot,u_beta)
load('FDx_Part.mat')  %@(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta)
load('FDy_Part.mat')  %@(Alpha,Alphadot,Alphaddot,Beta,Betadot,Betaddot,Theta,Thetadot,Thetaddot,u_beta)

p.FAx_Part = FAx_Part;         
p.FAy_Part = FAy_Part;
p.FDx_Part = FDx_Part;         
p.FDy_Part = FDy_Part;

load('X_Gene.mat')          %@(Alpha,Beta,Theta)
p.X_Gene = X_Gene;

load('f_beta_fn.mat')       %@(Alpha,Beta,Betadot,Theta)
load('g_beta_fn.mat')       %@(Alpha,Beta)

p.f_beta_fn = f_beta_fn;              p.g_beta_fn = g_beta_fn;        
end