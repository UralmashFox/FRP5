close all
clear all
clc
%syms A1 A2 A3
a1 = 1, a2 = 2, a3 = 3
t = 1
%- - - - - - - - - - - - - - - - - 
q1 = 0, q2 = 0, q3 = 0, dq1 = 1, ddq1 = 2, dq2 = 2, ddq2 = 1, dq3 = 0.5, ddq3 = 1.25
% - - - - - - - - - - - - - - - - -
m1 = 1, m2 = 2, m3 = 3, d1 = 8, d2 = 13
% syms q1 dq1 ddq1 q2 dq2 ddq2 q3 dq3 ddq3
%--------------------------------
while t<=350
    
    q1 = a1*sin(t)
    q2 = a2*cos(2*t)
    q3 = a3*sin(3*t)

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
 torque3 = - cross(f3, ([0 0 d2]+[0 0 q3/2])) + I3*domega3 + cross(omega3, I3*omega3)
 torque2 = torque3 - cross(f2, ([0 0 d1]+[0 0 d2/2])) + cross (f3, [0 d2/2 0]) + I2*domega2 + cross(domega2, I2*omega2)
 torque1 = torque2 - cross(f1, [0 0 d1/2]) + cross (f2, [0 0 d1/2]) + I1*domega1 + cross(domega1, I1*omega1)
 
 M3 = cross(f3, ([0 0 d2]+[0 0 q3/2]))
 M2 = torque3 - cross(f2, ([0 0 d1]+[0 0 d2/2]))
 M1 = torque2 - cross(f1, [0 0 d1/2])
 
 mo3(t) = sqrt(M3(1)^2+M3(2)^2+M3(3)^2)
 mo2(t) = sqrt(M2(1,1)^2+M2(2,1)^2+M2(3,1)^2 + M2(1,2)^2+M2(2,2)^2+M2(3,2)^2 + M2(1,3)^2+M2(2,3)^2+M2(3,3)^2)
 mo1(t) = sqrt(M1(1,1)^2+M1(2,1)^2+M1(3,1)^2 + M1(1,2)^2+M1(2,2)^2+M1(3,2)^2 + M1(1,3)^2+M1(2,3)^2+M1(3,3)^2)
 
 C3 = cross(omega3, I3*omega3)    
 C2 = cross(domega2, I2*omega2)
 C1 = cross(domega1, I1*omega1)
 
 c3(t) = sqrt(C3(1)^2+C3(2)^2+C3(3)^2)
 c2(t) = sqrt(C3(1)^2+C3(2)^2+C3(3)^2)
 c1(t) = sqrt(C3(1)^2+C3(2)^2+C3(3)^2)
 
 G3 = I3*domega3
 G2 = I2*domega2
 G1 = I1*domega1
 
 g3(t) = sqrt(G3(1)^2+G3(2)^2+G3(3)^2)
 g2(t) = sqrt(G3(1)^2+G3(2)^2+G3(3)^2)
 g1(t) = sqrt(G3(1)^2+G3(2)^2+G3(3)^2)
 
 tor3(t) = sqrt(torque3(1,1)^2+torque3(2,1)^2+torque3(3,1)^2 + torque3(1,2)^2+torque3(2,2)^2+torque3(3,2)^2 + torque3(1,3)^2+torque3(2,3)^2+torque3(3,3)^2)
 tor2(t) = sqrt(torque2(1,1)^2+torque2(2,1)^2+torque2(3,1)^2 + torque2(1,2)^2+torque2(2,2)^2+torque2(3,2)^2 + torque2(1,3)^2+torque2(2,3)^2+torque2(3,3)^2)
 tor1(t) = sqrt(torque1(1,1)^2+torque1(2,1)^2+torque1(3,1)^2 + torque1(1,2)^2+torque1(2,2)^2+torque1(3,2)^2 + torque1(1,3)^2+torque1(2,3)^2+torque1(3,3)^2)
 
 time(t) = t
 t = t+1
end
figure(1)
subplot(3,1,1)
plot(time,tor3)
title("\tau 3")
subplot(3,1,2)
plot(time,tor2)
title("\tau 2")
subplot(3,1,3)
plot(time,tor1)
title("\tau 1")

figure(2)
subplot(3,1,1)
plot(time,mo3)
title("M 3")
subplot(3,1,2)
plot(time,mo2)
title("M 2")
subplot(3,1,3)
plot(time,mo1)
title("M 1")

figure(3)
subplot(3,1,1)
plot(time,g3)
title("Gravity term 3")
subplot(3,1,2)
plot(time,g2)
title("Gravity term 2")
subplot(3,1,3)
plot(time,g1)
title("Gravity term 1")

figure(4)
subplot(3,1,1)
plot(time,c3)
title("Coriolis term 3")
subplot(3,1,2)
plot(time,c2)
title("Coriolis term 2")
subplot(3,1,3)
plot(time,c1)
title("Coriolis term 1")