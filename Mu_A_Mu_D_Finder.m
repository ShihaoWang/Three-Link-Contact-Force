function [mu_A,mu_D] = Mu_A_Mu_D_Finder(Post_Impact_Init,p,fai)

% This function is used to compute the offline database for the ICRA 2017
% conference
% This is the main function to compute the minimum friction coefficient

p.K = 10;   % This is the default gain that the joint motor will use

% Post-Impact integration
tspan = linspace(0,1.5,16);              % Assume that the system will be stabilize within 3s
options = odeset('AbsTol',1e-5,'RelTol',1e-5,'MassSingular','yes');

[t2,z2] = ode23(@rhs_post_impact_gain_K,tspan,Post_Impact_Init,options,p);

% However, if the initial values are at the boundary
% then the integration will go wrong.
if isnan(z2(end,1))==1
    
    mu_A = 10 + abs(rand);
    mu_D = 10 + abs(rand);
    
    % If the integration result is not trivial
else
    [mu_A,mu_D] = FriCoef_InteFor_Finder_BF(z2,p,fai);
end

end

