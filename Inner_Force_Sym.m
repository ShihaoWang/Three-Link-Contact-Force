function Inner_Force_Sym()
% This function is used to calculate the optimal inner force

syms A B C D E F G H M N x m n real

H_exp = sqrt((A + B * x)^2)/(C + D * x) + sqrt((E + F * x)^2)/(G + H * x);
% H_Neg = sqrt((A + B * x))/(C + D * x) - (E + F * x)/(G + H * x);

% Here we would like to take the derivate of this function with respect to
% the variable x to get the optimal value
% H = H_Pos;
% H = H_Neg;

H_Deri = jacobian(H_exp, x);
[H_Deri_N,~] = numden(H_Deri);  % Here what we got is a 2nd order polynomials y(x) = a x^2 + b x + c

H_Deri_N = simplify(H_Deri_N);

%% This next job is to make sure that the cubic equation can be solved in different cases
% collect(H_Deri_N, abs(A+B*x))
sgn_ApBX = - H*C^2*E^2 - H*C^2*E*F*x + G*C^2*E*F + G*C^2*F^2*x - 2*H*C*D*E^2*x - 2*H*C*D*E*F*x^2 + 2*G*C*D*E*F*x + 2*G*C*D*F^2*x^2 - H*D^2*E^2*x^2 - H*D^2*E*F*x^3 + G*D^2*E*F*x^2 + G*D^2*F^2*x^3;

% collect(H_Deri_N, abs(E+F*x))
sgn_EpFx = - D*A^2*G^2 - 2*D*A^2*G*H*x - D*A^2*H^2*x^2 - D*A*B*G^2*x + C*A*B*G^2 - 2*D*A*B*G*H*x^2 + 2*C*A*B*G*H*x - D*A*B*H^2*x^3 + C*A*B*H^2*x^2 + C*B^2*G^2*x + 2*C*B^2*G*H*x^2 + C*B^2*H^2*x^3;

eqn = sgn_ApBX * m + sgn_EpFx * n;

a = jacobian(jacobian(jacobian(eqn, x), x), x)/6;
b = (jacobian(jacobian(eqn, x), x) - 6 * a * x)/2;
c = jacobian(eqn, x) - 3 * a * x^2 - 2 * b * x;
d = eqn - a * x^3 - b * x^2 - c* x;
a = simplify(a);
b = simplify(b);
c = simplify(c);
d = simplify(d);

a_fn = matlabFunction(a);
b_fn = matlabFunction(b);
c_fn = matlabFunction(c);
d_fn = matlabFunction(d);

save('a_fn.mat', 'a_fn');
save('b_fn.mat', 'b_fn');
save('c_fn.mat', 'c_fn');
save('d_fn.mat', 'd_fn');

end

