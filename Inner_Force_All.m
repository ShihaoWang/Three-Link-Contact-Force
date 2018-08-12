function [x_tot, Feasi_tot, mu_A_tot, mu_D_tot] = Inner_Force_All(A, B, C, D, E, F, G, H)

x_low = ceil(max(-C/D, -G/H));
x_upp = x_low + 10;

x_tot = linspace(x_low, x_upp,101);

Feasi_tot = [];
mu_A_tot = [];
mu_D_tot = [];

for i = 1:length(x_tot)
    x = x_tot(i);
    [Feasi_i, mu_A_i, mu_D_i] = Fric_Coeff(A, B, C, D, E, F, G, H, x);
    Feasi_tot = [Feasi_tot; Feasi_i];
    mu_A_tot = [mu_A_tot; mu_A_i];
    mu_D_tot = [mu_D_tot; mu_D_i];
end
end

