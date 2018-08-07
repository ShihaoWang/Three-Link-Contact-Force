function Post_Impact_Cost_Database()

% This function is reserved for the computation of the offline database for
% the comprehensive case for ICRA 2017

p = Robot_Component_InertiaNLength();
p = Function_Attach(p);

load('thetadot_fn.mat')%@(Alpha,Beta,Betadot,Theta)
load('alphadot_fn.mat')%@(Alpha,Beta,Betadot,Theta)

N_Nor_Angl = 46;   % This is the number of points that the normal direction could take
N_Hz_Dist = 16;    % This is the number of points in the horizontal direction
N_Vt_Dist = 28;    % This is the number of points in the vertical direction
N_Beta = 41;       % This is the number of points for the beta angle
N_Betadot = 41;    % This is the number of points for the betadot

% Five main loops
Nor_Angl_array = linspace(0,pi/2, N_Nor_Angl);      % The polar angle of the contact point is from 0 to pi/2
Hz_Dist_array = linspace(0.1, 0.25, N_Hz_Dist);     % The horizontal distance is from 10cm to 25cm
Vt_Dist_array = linspace(0, 0.27,N_Vt_Dist);        % The vertical distance is 0cm to 28cm
Beta_array = linspace(pi/4, pi,N_Beta);             % Here Beta stands for the Post-Impact Beta
Betadot_array = linspace(-6,6,N_Betadot);           % Here Betadot stands for the Post-Impact Betadot

Nor_Angl_unit =  Nor_Angl_array(2) - Nor_Angl_array(1);
Hz_Dist_unit =   Hz_Dist_array(2) - Hz_Dist_array(1);
Vt_Dist_unit =   Vt_Dist_array(2) - Vt_Dist_array(1);
Beta_unit =      Beta_array(2) - Beta_array(1);
Betadot_unit =   Betadot_array(2) - Betadot_array(1);

% Mu_A/E_lib saves the friction coefficient for contact A and contact E
Mu_A_lib = zeros(N_Nor_Angl, N_Hz_Dist,N_Vt_Dist, N_Beta,N_Betadot);
Mu_D_lib = zeros(N_Nor_Angl, N_Hz_Dist,N_Vt_Dist, N_Beta,N_Betadot);

% Save em into 1-D library
One_D_Mu_A_lib = [];
One_D_Mu_D_lib = [];

for i = 46:N_Nor_Angl    % Iteration for the normal angle of the contact point
    
    Nor_Angl_i = Nor_Angl_array(1) + (i-1) * Nor_Angl_unit;
    
    for j = 1:N_Hz_Dist % Iteration for different  horizontal distances
        
        Hz_Dist_j = Hz_Dist_array(1) + (j-1) * Hz_Dist_unit;
        
        for k = 1: N_Vt_Dist % Iteration for different vertical distances
            
            Vt_Dist_k = Vt_Dist_array(1) + (k - 1) * Vt_Dist_unit;
            
            for l = 1:N_Beta % Iteraction for different Betas
                
                Beta_l = Beta_array(1) + (l - 1) * Beta_unit;
                
                for m = 1:N_Betadot  % Iteraction for different Betadots
                    
                    Betadot_m = Betadot_array(1) + (m - 1) *Betadot_unit;
                    
                    %-----------------------------------------------------%
                    fai = Nor_Angl_i;
                    Wall_x = Hz_Dist_j;
                    Wall_y = Vt_Dist_k;
                    Beta = Beta_l;
                    Betadot = Betadot_m;
                    
                    % Post-Impact anlge: theta
                    theta1 = 2*atan((10000*Wall_y + (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 10000*Wall_x + 40000*Wall_y^2 + 560*cos(Beta) + 29));
                    theta2 = 2*atan((10000*Wall_y - (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 10000*Wall_x + 40000*Wall_y^2 + 560*cos(Beta) + 29));
                    
                    % Post-Impact angle: alpha
                    alpha1 = -2*atan((1000*sin(Beta) + (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 40000*Wall_y^2 + 1560*cos(Beta) - 1921));
                    alpha2 = -2*atan((1000*sin(Beta) - (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 40000*Wall_y^2 + 1560*cos(Beta) - 1921));
                    
                    % Nontrivial solution will form into a closure chain
                    if (isreal(theta1))&&(isreal(alpha1))                     
                        % Choose the bigger Theta
                        Theta = max(theta1, theta2);
                        if (Theta == theta1)
                            Alpha = alpha1;
                        else
                            Alpha = alpha2;
                        end
                        
                        Thetadot = thetadot_fn(Alpha,Beta,Betadot);
                        Alphadot = alphadot_fn(Alpha,Beta,Betadot);
                        Post_Impact_Init = [Theta;...
                            Alpha;...
                            Beta;...
                            Thetadot;...
                            Alphadot;...
                            Betadot];
                        
                        [mu_A,mu_D] = Mu_A_Mu_D_Finder(Post_Impact_Init,p,fai);
                        
                        Mu_A_lib(i,j,k,l,m) = mu_A;
                        mu_D_lib(i,j,k,l,m) = mu_D;
                        
                        One_D_Mu_A_lib = [One_D_Mu_A_lib;mu_A];
                        One_D_Mu_D_lib = [One_D_Mu_D_lib;mu_D];
                    else                      
                        Mu_A_lib(i,j,k,l,m) = 10 + abs(rand);
                        Mu_D_lib(i,j,k,l,m) = 10 + abs(rand);
                        
                        One_D_Mu_A_lib = [One_D_Mu_A_lib;mu_A];
                        One_D_Mu_D_lib = [One_D_Mu_D_lib;mu_D];
                    end
                    %-----------------------------------------------------%
                end
            end
        end
    end
end

end

