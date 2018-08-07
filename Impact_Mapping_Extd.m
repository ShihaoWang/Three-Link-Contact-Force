function [Mag,Post_Impact_State]= Impact_Mapping_Extd(z)

Theta = z(1);     Alpha = z(2);     Beta = z(3); 
Thetadot = z(4);  Alphadot = z(5);  Betadot = z(6);

load('Pre_T_fn');

Pre_T_Val = Pre_T(Alpha,Alphadot,Beta,Betadot,Thetadot);
Post_Impact_Betadot = sign(Betadot) * sqrt(Pre_T_Val/(2*(319216919008673693160713*cos(2*Alpha) + 356978521001049373636224*cos(2*Beta) + 6089075567150829707155456*cos(2*Alpha + 2*Beta) - 6765271007160552773952393)/(4611686018427387904000000*(49*cos(2*Alpha) + 100*cos(2*Alpha + 2*Beta) + 140*cos(Beta) - 140*cos(2*Alpha + Beta) - 149))));
Post_Impact_Thetadot = -(28*Betadot*sin(Beta))/(5*(10*sin(Alpha + Beta) - 7*sin(Alpha)));
Post_Impact_Alphadot = -(2*Betadot*(25*sin(Alpha + Beta) - 14*sin(Beta)))/(5*(10*sin(Alpha + Beta) - 7*sin(Alpha)));


Linear_Mom_Coef = [ (5523*sin(Alpha + Theta))/250000 - sin(Alpha + Beta + Theta)/400 - (12827*sin(Theta))/160000, (5523*sin(Alpha + Theta))/250000 - sin(Alpha + Beta + Theta)/400, -sin(Alpha + Beta + Theta)/400;
    cos(Alpha + Beta + Theta)/400 - (5523*cos(Alpha + Theta))/250000 + (12827*cos(Theta))/160000, cos(Alpha + Beta + Theta)/400 - (5523*cos(Alpha + Theta))/250000,  cos(Alpha + Beta + Theta)/400];

Mag = norm(Linear_Mom_Coef*([Post_Impact_Thetadot;Post_Impact_Alphadot;Post_Impact_Betadot] - [Thetadot; Alphadot; Betadot]));
% load('Impulse_De_fn.mat');   %@(Alpha,Beta,Theta)
% load('Impulse_E_fn.mat');   %@(Alpha,Beta,Theta)
% 
% De = Impulse_De_fn(Alpha,Beta,Theta);
% E = Impulse_E_fn(Alpha,Beta,Theta);
% 
% f = -(E * De^(-1) * E')^(-1) * E * [Thetadot; Alphadot; Betadot; 0;0];
% 
% Impulse_Ax = f(1);
% Impulse_Ay = f(2);
% Impulse_Ex = f(3);
% Impulse_Ey = f(4);
% 
% Impulse_x = Impulse_Ax + Impulse_Ex;
% Impulse_y = Impulse_Ay + Impulse_Ey;
% 
% Mag = sqrt(Impulse_x * Impulse_x + Impulse_y * Impulse_y);
% 
% post_int_extd = [Thetadot; Alphadot; Betadot; 0; 0] + De^(-1)*(E' *f);
% 
% Post_Impact_State = zeros(6,1);

Post_Impact_State(1) = Theta;
Post_Impact_State(2) = Alpha;
Post_Impact_State(3) = Beta;
Post_Impact_State(4) = Post_Impact_Thetadot;
Post_Impact_State(5) = Post_Impact_Alphadot;
Post_Impact_State(6) = Post_Impact_Betadot;

end