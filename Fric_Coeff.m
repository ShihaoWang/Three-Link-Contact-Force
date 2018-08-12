function [Feasi, mu_A, mu_D] = Fric_Coeff(A, B, C, D, E, F, G, H, x)

A_Num = A + B * x;
A_Den = C + D * x;

D_Num = E + F * x;
D_Den = G + H * x;

if (A_Den>0)&&(D_Den>0)
    Feasi = 1;
else
    Feasi = 0;
end

if Feasi == 1
    mu_A = abs(A_Num/A_Den);
    mu_D = abs(D_Num/D_Den);
else
    mu_A = 1;
    mu_D = 1;
end

end