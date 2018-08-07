function InteFor_Analytic_Comp()

% This function is used to compute the analytic expression for the
% corresponding c that can bring the sum of the two coefficients to achieve
% the minimum value

syms a b c d e f g h x fai real

mu_A = ( a + x * b )/( c + x * d );
% mu_E = ( e - x * (b * sin(fai) + d * cos(fai)) )/( f + x * (d * sin(fai) - cos(fai) * b));

mu_E = (e + x * f)/(g + x * h);

mu_fn = abs(mu_A) + abs(mu_E);

Jac_x = simplify(jacobian(mu_fn,x));

Jac_x_extd = simplifyFraction(Jac_x * (c + d*x)^2 * (g + h * x)^2);

% collect(Jac_x_extd);
% (2*c*b^2*h^3 - 2*a*b*d*h^3 + 2*g*d^3*f^2 - 2*e*d^3*f*h)*x^4 + 
% (- 2*a^2*d*h^3 - 6*g*a*b*d*h^2 + 2*c*a*b*h^3 + 6*c*g*b^2*h^2 - 2*d^3*e^2*h + 2*g*d^3*e*f - 6*c*d^2*e*f*h + 6*c*g*d^2*f^2)*x^3 + 
% (- 6*a^2*d*g*h^2 + 6*a*b*c*g*h^2 - 6*a*b*d*g^2*h + 6*b^2*c*g^2*h - 6*c^2*d*e*f*h + 6*c^2*d*f^2*g - 6*c*d^2*e^2*h + 6*c*d^2*e*f*g)*x^2 + 
% (- 6*d*h*a^2*g^2 + 6*h*a*b*c*g^2 - 2*d*a*b*g^3 + 2*b^2*c*g^3 - 2*h*c^3*e*f + 2*c^3*f^2*g - 6*d*h*c^2*e^2 + 6*d*c^2*e*f*g)*x 
% - 2*d*a^2*g^3 + 2*b*a*c*g^3 - 2*h*c^3*e^2 + 2*f*c^3*e*g

x_sym = solve(Jac_x,x,'MaxDegree',4);

x_sym = simplify(x_sym)
end

