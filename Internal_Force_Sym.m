function Internal_Force_Sym()
% This function is used to calculate the optimal inner force

syms A B C D E F G H M N x real

H_x = ((A + B * x)/(C + D * x))^2 + ((E + F * x)/(G + H * x))^2;

H_Deri = jacobian(H_x, x);

[H_Deri_N,~] = numden(H_Deri);  % Here what we got is a 2nd order polynomials y(x) = a x^2 + b x + c
H_Deri_N = simplify(H_Deri_N);

% H_Deri_N_fn = matlabFunction(H_Deri_N);
% save('H_Deri_N_fn.mat','H_Deri_N_fn');


X = fzero(@(x) H_Deri_N, 0) 



% a = simplify(jacobian(jacobian(H_Deri_N, x),x)/2);
% b = simplify(jacobian(H_Deri_N,x) - 2 * a * x);
% c = simplify(H_Deri_N - a * x^2 - b * x);
% 
% r = simplify(roots([ a b c ]));
% 
% % r_Pos_fn = matlabFunction(r);
% r_Neg_fn = matlabFunction(r);
% 
% % save('r_Pos_fn.mat','r_Pos_fn')
% save('r_Neg_fn.mat','r_Neg_fn')

end

