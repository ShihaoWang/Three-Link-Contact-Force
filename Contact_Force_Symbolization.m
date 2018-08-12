function [A, B, C, D, E, F, G, H] = Contact_Force_Symbolization(Part_Tot, Contact_X, Contact_Y, Envi_Theta)

FAx_Par = Part_Tot(1);          FAy_Par = Part_Tot(2);              FDx_Par = Part_Tot(3);          FDy_Par = Part_Tot(4);

Inner_Direction = [Contact_X Contact_Y]';                           Inner_Direction = Inner_Direction/norm(Inner_Direction);

A = FAx_Par;                    C = FAy_Par;                        

FD_Par = [FDx_Par; FDy_Par];

[FDn_Par, FDt_Par] = Force_Proj(FD_Par, Envi_Theta);

[FDn_Gen, FDt_Gen] = Force_Proj(-Inner_Direction, Envi_Theta);

B = Inner_Direction(1);

D = Inner_Direction(2);

E = FDt_Par;            F = FDt_Gen;
G = FDn_Par;            H = FDn_Gen;

end

function [Force_Norm, Force_Tang] = Force_Proj(Force, Theta)

Force_Norm = Force(2) * cos(Theta) - Force(1) * sin(Theta);                                                 
Force_Tang = Force(1) * cos(Theta) + Force(2) * sin(Theta);

end

