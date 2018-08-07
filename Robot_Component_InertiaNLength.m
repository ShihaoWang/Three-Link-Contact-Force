function p = Robot_Component_InertiaNLength()

% This function is used to compute the corresponding moment of inertia from
% each component of the robot to the three-link model

m_body = 345.9/1000;  
L_body = 10.5/100;      W_body = 7/100;        H_body = 6/100;  
I_body = MOI(m_body, L_body, W_body, H_body);    

m_shoulder = 7.5/1000;  
L_shoulder = 6/100;     W_shoulder = 2.5/100;  H_shoulder = 3.5/100;  
I_shoulder = MOI(m_shoulder, L_shoulder, W_shoulder, H_shoulder);  

m_armNhand = 50/1000;
L_armNhand = 10/100;    W_armNhand = 2.5/100;  H_armNhand = 3.5/100;
I_armNhand = MOI(m_armNhand, L_armNhand, W_armNhand, H_armNhand);  

m_pelvis = 22.1/1000;
L_pelvis = 3.5/100;     W_pelvis = 3/100;      H_pelvis = 3/100;
I_pelvis = MOI(m_pelvis, L_pelvis, W_pelvis, H_pelvis);  

m_leg = 55.6/1000;
L_leg = 12/100;         W_leg = 3.5/100;        H_leg = 3/100;
I_leg = MOI(m_leg, L_leg, W_leg, H_leg);  

m_foot = 42.6/1000;
L_foot = 1/100;         W_foot = 4.5/100;       H_foot = 8.5/100;
I_foot = MOI(m_foot, L_foot, W_foot, H_foot);  

% First link is from the foot edge to the pelvis
p.m1 = m_foot + m_leg + m_pelvis;
p.L1 = 12.5/100;
p.I1 = I_leg + (I_pelvis + m_pelvis * L_leg * L_leg/4) + (I_foot + m_foot * L_leg * L_leg/4);

% Second link is from the pelvis to the shoulder/arm joint
p.m2 = m_body + m_pelvis + m_leg + m_foot + 2*m_shoulder + m_armNhand;
p.L2 = 7/100;
p.I2 = (I_body + m_body * (W_body^2/4))    + (I_pelvis + m_pelvis * (0.06^2 + 0.05^2)) +...
       (I_leg + m_leg * (0.65^2 + 0.09^2)) + (I_foot + m_foot * (0.085^2 + 0.15^2)) + ...
       (2*I_shoulder + m_shoulder * (0.1^2+ 0.03^2 + 0.03^2 + 0.01^2)) + (I_armNhand + m_armNhand * (0.11^2 + 0.01^2));

% Third link is from the shoulder/arm joint to the hand tip
p.m3 = m_armNhand;
p.L3 = 10/100;
p.I3 = I_armNhand;


% The maximum angular velocity is 7.3 rad/s
% The maximum angular acceleration is 104 rad/s^2
% The ratio between these two are 104/7.4 = 14
% p.K = 14;
end


function moment_of_inertia = MOI(m,L,W,H)
moment_of_inertia = 1/12 * m * ( L^2 + W^2);
end

