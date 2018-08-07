function [mu_A,mu_D] = FriCoef_InteFor_Finder_BF(Config_Data,p,fai)

% This function will compute the the current minimum mu_A and mu_D

% Load the generalize force
load('X_Gene');

% Particular solutions
load('FAx_Part.mat');
load('FAy_Part.mat');
load('FDx_Part.mat');
load('FDy_Part.mat');

% The Post-Impact State
Theta_array = Config_Data(:,1);      Alpha_array = Config_Data(:,2);       Beta_array = Config_Data(:,3);
Thetadot_array = Config_Data(:,4);   Alphadot_array = Config_Data(:,5);    Betadot_array = Config_Data(:,6);

if (sin(Alpha_array(1) + Beta_array(1) + Theta_array(1))/10 - (7*sin(Alpha_array(1) + Theta_array(1)))/100 + sin(Theta_array(1))/8) < 0.01    
    
    FAx_Gene = 1; 
    FAy_Gene = 0;
    FDx_Gene = -1;
    FDy_Gene = 0;
    
else
    % %Generalised force
    Gene_fn = -X_Gene(Alpha_array(1),Beta_array(1),Theta_array(1));
    
    % General forces
    FAx_Gene = Gene_fn(1);
    FAy_Gene = Gene_fn(2);
    FDx_Gene = Gene_fn(3);
    FDy_Gene = Gene_fn(4);
    
    % Sometimes it is needed to normalize the vector to make it of magnitude 1
    FAx_Gene = FAx_Gene/sqrt(Gene_fn(1)^2 + Gene_fn(2)^2);
    FAy_Gene = FAy_Gene/sqrt(Gene_fn(1)^2 + Gene_fn(2)^2);
    FDx_Gene = FDx_Gene/sqrt(Gene_fn(3)^2 + Gene_fn(4)^2);
    FDy_Gene = FDy_Gene/sqrt(Gene_fn(3)^2 + Gene_fn(4)^2);
end

% Generalized force projection
FDN_Gene = cos(fai) * FDx_Gene - sin(fai) * FDy_Gene;  % This force pointing inward the normal contact plane
FDT_Gene = sin(fai) * FDx_Gene + cos(fai) * FDy_Gene;  % Thie force is in the normal contact plane 

Force_Gene = [FAx_Gene; FAy_Gene;FDN_Gene; FDT_Gene];

A_fn = p.Post_A_fn;
B_fn = p.Post_B_fn;

g1_fn = p.g1_fn;
f1_fn = p.f1_fn;

K = p.K;

% Initialize the friction coefficients
mu_A_array = [];
mu_D_array = [];

for i = 1:length(Theta_array)
    
    % The state for each time stpe
    Theta_i = Theta_array(i);          Alpha_i = Alpha_array(i);            Beta_i = Beta_array(i);
    
    Thetadot_i = Thetadot_array(i);    Alphadot_i = Alphadot_array(i);      Betadot_i = Betadot_array(i);
    
    g1_fn_i = g1_fn(Alpha_i,Beta_i);
    
    f1_fn_i = f1_fn(Alpha_i,Beta_i,Betadot_i,Theta_i);
    
    % The control torque for each time step
    u_beta_i = (-K * Betadot_i - f1_fn_i)/g1_fn_i;
    
    A_i = A_fn(Alpha_i,Beta_i,Theta_i);
    
    B_i = B_fn(Alpha_i,Alphadot_i,Beta_i,Betadot_i,Theta_i,Thetadot_i,u_beta_i);
    
    stateddot_i = -A_i\B_i;
    % Stateddot
    Thetaddot_i = stateddot_i(1);    Alphaddot_i = stateddot_i(2);    Betaddot_i = stateddot_i(3);
    
    % Post-Impact particular forces
    FAx_Part_i = FAx_Part(Alpha_i,Alphadot_i,Alphaddot_i,Beta_i,Betaddot_i,Theta_i,Thetadot_i,Thetaddot_i,u_beta_i);
    FAy_Part_i = FAy_Part(Alpha_i,Alphadot_i,Alphaddot_i,Beta_i,Betaddot_i,Theta_i,Thetadot_i,Thetaddot_i,u_beta_i);
    FDx_Part_i = FDx_Part(Alpha_i,Alphadot_i,Alphaddot_i,Beta_i,Betadot_i,Betaddot_i,Theta_i,Thetadot_i,Thetaddot_i,u_beta_i);
    FDy_Part_i = FDy_Part(Alpha_i,Alphadot_i,Alphaddot_i,Beta_i,Betadot_i,Betaddot_i,Theta_i,Thetadot_i,Thetaddot_i,u_beta_i);
    
    % Particular Force projection
    FDN_Part_array_i = cos(fai) * FDx_Part_i - sin(fai) * FDy_Part_i;
    FDT_Part_array_i = sin(fai) * FDx_Part_i + cos(fai) * FDy_Part_i;
    
    Force_Part = [FAx_Part_i;FAy_Part_i;FDN_Part_array_i; FDT_Part_array_i];
    
    [mu_A_i, mu_D_i] = FriCoef_Cur_i(Force_Part,Force_Gene);
  
    mu_A_array = [mu_A_array;mu_A_i];
    mu_D_array = [mu_D_array;mu_D_i];

end

mu_sum = mu_A_array + mu_D_array;

[mu_max, mu_ind] = max(mu_sum);

mu_A = mu_A_array(mu_ind);

mu_D = mu_D_array(mu_ind);

end

function [mu_A, mu_D] = FriCoef_Cur_i(Force_Part,Force_Gene)

% This funciton computes the mu_A and mu_D given the current state

% Particular Forces
FAx_p = Force_Part(1);
FAy_p = Force_Part(2);
FDN_p = Force_Part(3);
FDT_p = Force_Part(4);

% Gneralize Forces
FAx_g = Force_Gene(1);
FAy_g = Force_Gene(2);
FDN_g = Force_Gene(3);
FDT_g = Force_Gene(4);

Force_Part = [FAx_p; FAy_p; FDN_p;FDT_p];
Force_Gene = [FAx_g; FAy_g; FDN_g;FDT_g];

[mu_A, mu_D]= Mu_Selec(Force_Part, Force_Gene);

end

function c_low = C_Low_Bound(Force_Part,Force_Gene)

% This function calculates the low bound of c
% Particular Forces
FAy_p = Force_Part(2);
% Gneralize Forces
FAy_g = Force_Gene(2);

c_low = -FAy_p/FAy_g;
c_low = ceil(c_low);

end

function [mu_A, mu_D]= Mu_Selec(Force_Part,Force_Gene)
% The low bound to gurantee the validity of the c* value
c_low = C_Low_Bound(Force_Part,Force_Gene);

a = Force_Part(1);     % FAx_p  ==> a
c = Force_Part(2);     % FAy_p  ==> c
g = Force_Part(3);     % FDN_p  ==> g
e = Force_Part(4);     % FDT_p  ==> e

b = Force_Gene(1);     % FAx_g  ==> b
d = Force_Gene(2);     % FAy_g  ==> d
h = Force_Gene(3);     % FDN_g  ==> h
f = Force_Gene(4);     % FDT_g  ==> f

p4 = (2*c*b^2*h^3 - 2*a*b*d*h^3 + 2*g*d^3*f^2 - 2*e*d^3*f*h);
p3 = (- 2*a^2*d*h^3 - 6*g*a*b*d*h^2 + 2*c*a*b*h^3 + 6*c*g*b^2*h^2 - 2*d^3*e^2*h + 2*g*d^3*e*f - 6*c*d^2*e*f*h + 6*c*g*d^2*f^2);
p2 = (- 6*a^2*d*g*h^2 + 6*a*b*c*g*h^2 - 6*a*b*d*g^2*h + 6*b^2*c*g^2*h - 6*c^2*d*e*f*h + 6*c^2*d*f^2*g - 6*c*d^2*e^2*h + 6*c*d^2*e*f*g);
p1 = (- 6*d*h*a^2*g^2 + 6*h*a*b*c*g^2 - 2*d*a*b*g^3 + 2*b^2*c*g^3 - 2*h*c^3*e*f + 2*c^3*f^2*g - 6*d*h*c^2*e^2 + 6*d*c^2*e*f*g);
p0 = - 2*d*a^2*g^3 + 2*b*a*c*g^3 - 2*h*c^3*e^2 + 2*f*c^3*e*g;

p = [p4 p3 p2 p1 p0];

x_opt_array = roots(p);   % From the polynomial, we know that there are four solutions 
                          % But only two solutions are nontrivial                     
mu_sum_array = 10 * ones(4,1);

for i = 1:length(x_opt_array)
    x_opt_i = x_opt_array(i);
    if isreal(x_opt_i)
        if x_opt_i>=c_low
            mu_A_i = abs((a + x_opt_i * b)/(c + x_opt_i * d));
            mu_D_i = abs((e + x_opt_i * f)/(g + x_opt_i * h));
            mu_sum_array(i) = mu_A_i + mu_D_i;           
        end
    end
end

[mu_min,mu_ind] = min(mu_sum_array);

if mu_min<10    
    x_opt = x_opt_array(mu_ind);
    mu_A = abs((a + x_opt * b)/(c + x_opt * d));
    mu_D = abs((e + x_opt * f)/(g + x_opt * h)); 
    
else
    mu_A = 10 + abs(rand);
    mu_D = 10 + abs(rand);
    
end

end