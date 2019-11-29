close all
clear all
clc
%- - - - - - - - - - - - - - - - - 
q1 = deg2rad(90), dq1 = 1, ddq1 = 2, q2 = deg2rad(90)
dq2 = 2, ddq2 = 1, q3 = 5, dq3 = 0.5, ddq3 = 1.25
% - - - - - - - - - - - - - - - - -
m1 = 1, m2 = 2, m3 = 3, d1 = 8, d2 = 13
% syms q1 dq1 ddq1 q2 dq2 ddq2 q3 dq3 ddq3
%--------------------------------
R0  = [cos(q1) -sin(q1) 0
       sin(q1)  cos(q1) 0
           0       0    1]
    
Rt0 = [1    0   0
       0    0  -1
       0    1   0]
   
 A0 = R0*Rt0
%----------------------------------
 R1  = [cos(q2) -sin(q2) 0
        sin(q2)  cos(q2) 0
            0       0    1]
 
 Rt1 = [0   0   1
        1   0   0
        0   1   0]
    
 A1 = R1*Rt1  
 %---------------------------------
 Rt2 = [1    0   0
        0    1   0
        0    0   1]
    
 A2 = Rt2
 %--------Angular velocity-----------------
 % First (revolute) wrt base
 omega1 = [1 0 0; 0 1 0; 0 0 1]' * (ddq1*[0 0 1])'
 % Second (revolute) wrt base
 omega2 = A0'*(omega1 + (ddq2*[0 0 1])')
 % Third (prismatic) wrt base
 omega3 = A1*(omega1 + (ddq2*[0 0 1])')
 OMEGA = [omega1 omega2 omega3]
 %---------Angular acceleration-------------
 % First (revolute) wrt base
 domega1 = [1 0 0; 0 1 0; 0 0 1]*(ddq1*[0 0 1]' + cross(dq1*[0 0 0], [0 0 1])')
 % Second (revolute) wrt base
 domega2 = A0'*(domega1 + ddq2*[0 0 1]' + cross(dq2*(ddq1*[0 0 1]' + cross(dq1*[0 0 0], [0 0 1])'), [0 0 1])') 
 % Third (prismatic) wrt base
 domega3 = A1'*domega2 
 dOMEGA = [domega1 domega2 domega3]
 %---------Linear velocity-------------------
 % First (revolute) wrt base
 v1 = 0 + cross(omega1, [0 0 0])
 % Second (revolute) wrt base
 v2 = v1 + cross(omega2, [0 0 d1]) 
 % Third (prismatic) wrt base
 v3 = v2 + dq3*[0 0 1] + cross(omega3, [0 0 d2])
 V = [v1 v2 v3]
 %---------Linear acceleration-------------------
 % First (revolute) wrt base
 dv1 = cross(domega1, [0 0 0])
 % Second (revolute) wrt base
 dv2 = dv1 + cross(domega2, [0 0 d1]) + cross(omega2, cross(omega2, [0 0 d1]))
 % Third (prismatic) wrt base
 dv3 = dv2 + ddq3*[0 0 1] + cross(2*dq3*omega3, [0 0 1]) + cross(domega1, [0 0 d2]) + cross(omega3, cross(omega3, [0 0 d2]))
 dV = [dv1 dv2 dv3]
 %Acceleration of center of mass:
 dvc1 = dv1 + cross(domega1, [0 0 1/2*d1]) + cross(omega1, cross(omega1, [0 0 1/2*d1]))
 dvc2 = dv2 + cross(domega2, [0 0 1/2*d2]) + cross(omega2, cross(omega2, [0 0 1/2*d2]))
 dvc3 = dv3 + cross(domega3, [0 0 1/2*q3]) + cross(omega3, cross(omega3, [0 0 1/2*q3]))
 %Backward recursion: forces
 f3 = m3*dvc3
 f2 = f3 + m2*dvc2
 f1 = f2 + m1*dvc1
 %Moment of inertia:
 I1 = 1/12*m1*d1^2
 I2 = 1/12*m2*d2^2
 I3 = 1/12*m3*q3^2
 %Backward recursion: torques
 t3 = - cross(f3, ([0 0 d2]+[0 0 q3/2])) + I3*domega3 + cross(omega3, I3*omega3)
 t2 = t3 - cross(f2, ([0 0 d1]+[0 0 d2/2])) + cross (f3, [0 d2/2 0]) + I2*domega2 + cross(domega2, I2*omega2)
 t1 = t2 - cross(f1, [0 0 d1/2]) + cross (f2, [0 0 d1/2]) + I1*domega1 + cross(domega1, I1*omega1)
 
 
 
 