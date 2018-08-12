function [ mu_A, mu_D, x] = Inner_Force_Optimization(A, B, C, D, E, F, G, H, p)

% This function is used to calculate the optimal inner force magnitude

options = optimoptions('fmincon','Display','off','Algorithm','sqp');

p.A = A;            
p.B = B;            
p.C = C;            
p.D = D; 
p.E = E;            
p.F = F;            
p.G = G;            
p.H = H; 

[x,fval] = fmincon(@Internal_Force_Opt,0,[],[],[],[],[],[],@Internal_Force_Con, options, p);

mu_A = abs(A + B * x)/(C + D * x);
mu_D = abs(E + F * x)/(G + H *x);
end

function Obj = Internal_Force_Opt(x, p)

A = p.A;            B = p.B;            C = p.C;            D = p.D; 
E = p.E;            F = p.F;            G = p.G;            H = p.H; 

Obj = abs((A + B * x)/(C + D * x)) + abs((E + F * x)/(G + H *x));

end


function [c,ceq] = Internal_Force_Con(x,p)

eps_den = 0.1;
A = p.A;            B = p.B;            C = p.C;            D = p.D; 
E = p.E;            F = p.F;            G = p.G;            H = p.H; 
c = [-(C + D * x) + eps_den];
c = [c; -(G + H * x) + eps_den];

% c = [c ; -abs(A + B * x)/(C + D * x)];
% c = [c ; -1 + abs(A + B * x)/(C + D * x)];
% c = [c ; -abs(E + F * x)/(G + H *x)];
% c = [c ; -1 + abs(E + F * x)/(G + H *x)];

ceq = [];

end

