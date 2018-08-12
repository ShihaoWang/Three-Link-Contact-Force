function [mu_A, mu_D, x] = Inner_Force_Analytic(A, B, C, D, E, F, G, H, p)

% This function is used to calculate the analytic optimal solution for the  friction coefficient problem
a_fn = p.a_fn;          b_fn = p.b_fn;          c_fn = p.c_fn;          d_fn = p.d_fn;

% There are four possible cases:
% case 1: m = 1, n = 1
m = 1;      n = 1;
a = a_fn(A,B,C,D,E,F,G,H,m,n);
b = b_fn(A,B,C,D,E,F,G,H,m,n);
c = c_fn(A,B,C,D,E,F,G,H,m,n);
d = d_fn(A,B,C,D,E,F,G,H,m,n);
x = roots([a b c d]);

% case 2: m = 1, n = -1
m = 1;      n = -1;
a = a_fn(A,B,C,D,E,F,G,H,m,n);
b = b_fn(A,B,C,D,E,F,G,H,m,n);
c = c_fn(A,B,C,D,E,F,G,H,m,n);
d = d_fn(A,B,C,D,E,F,G,H,m,n);
x = [x; roots([a b c d])];

% case 3: m = -1, n = 1
m = -1;      n = 1;
a = a_fn(A,B,C,D,E,F,G,H,m,n);
b = b_fn(A,B,C,D,E,F,G,H,m,n);
c = c_fn(A,B,C,D,E,F,G,H,m,n);
d = d_fn(A,B,C,D,E,F,G,H,m,n);
x = [x; roots([a b c d])];

% case 4: m = -1, n = -1
m = -1;      n = -1;
a = a_fn(A,B,C,D,E,F,G,H,m,n);
b = b_fn(A,B,C,D,E,F,G,H,m,n);
c = c_fn(A,B,C,D,E,F,G,H,m,n);
d = d_fn(A,B,C,D,E,F,G,H,m,n);
x = [x; roots([a b c d])];

mu_A_tot = [];
mu_D_tot = [];

% Now it is time to characterize these solution
x_tot = [];

for i = 1:length(x)
    x_i = x(i);
    % First test is the real/imaginary test
    if isreal(x_i) == 1
        x_tot = [x_tot; x_i];
        % Second test is the normal force feasibility test
        Norm_A = C + x_i * D;
        Norm_D = G + x_i * H;
        if (Norm_A>0)&&(Norm_D>0)
            mu_A_i = abs((A + B * x_i)/Norm_A);
            mu_D_i = abs((E + F * x_i)/Norm_D);
            mu_A_tot = [mu_A_tot; mu_A_i];
            mu_D_tot = [mu_D_tot; mu_D_i];
        end
        
    end
end

mu_sum_tot = mu_A_tot + mu_D_tot;
[~, Min_Ind] = min(mu_sum_tot);
x = x_tot(Min_Ind);
mu_A = mu_A_tot(Min_Ind);
mu_D = mu_D_tot(Min_Ind);

end