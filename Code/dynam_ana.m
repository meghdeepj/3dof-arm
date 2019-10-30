%function for finding joint parameters for given trajectory of motion
function [T1 T2 T3 q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 t] = dynam_ana
 
%mass, moment of intertia, lengths and centre of mass distance for links
m1 = 0.14; m2=0.21; m3=0.07;
I1 = 4.7133e-04; I2=0.0016; I3=6.0667e-05;
l1 = 0.2; l2=0.3; l3=0.1;
r1= 0.1; r2 =0.15; r3=0.05;

%time taken from 0 to 4s for 1 revolution
%time step of 0.01secs
t= 0:0.01:4;
t=t';

%trajectory equation 
%origin taken as joint 1
x = 0.4+0.04*cos(pi*t/2);
y = 0.1+0.04*sin(pi*t/2);

%setting up parameters for fsolve function for finding joint parameters for
%joint parameters for each time step
x0 = [0 0];
temp = zeros(2,401);  

%--------------------------------------------------------------------------
%loop for finding joint angles in time steps of 0.01 from trajectory eqn
%relating x and y co-ordinates with joint angles
for i = 1:401   
    a = x(i);
    b = y(i);
        f = @(q)[0.2*cos(q(1))+0.3*cos(q(1)+q(2))+0.1-a , 0.2*sin(q(1))+0.3*sin(q(1)+q(2))-b];
        [temp(:,i)] = fsolve(f,x0);
end

%--------------------------------------------------------------------------
%extracting join angles from temp matrix
q1 = temp(1,:)';
q2 = 2*pi+temp(2,:)';
q3 = 2*pi-(q1+q2);

%finding joint parameter - angular velocity in time steps
dq1 = diff(q1)./diff(t);
dq1 = [dq1(1); dq1];
dq2 = diff(q2)./diff(t);
dq2 = [dq2(1); dq2];
dq3 = diff(q3)./diff(t);
dq3 = [dq3(1); dq3];

%finding joint parameter - angular acceleration in time steps
ddq1 = diff(dq1)./diff(t);
ddq1 = [ddq1(1); ddq1];
ddq2 = diff(dq2)./diff(t);
ddq2 = [ddq2(1); ddq2];
ddq3 = diff(dq3)./diff(t);
ddq3 = [ddq3(1); ddq3];

%--------------------------------------------------------------------------
%using derived equations of motion to find joint parameter - torques in
%torques at each joint time steps
T1 = (m1*r1^2 + 2*m2*l1*r2*cos(q1)+ m2*r2^2+ m3*l1^2 + 2*m3*l1*l2*cos(q2) + 2*m3*l1*r3*cos(q1+q2) + m3*l2^2 + 2*m3*l2*r3*cos(q3) + I1 + I2 + I3) .* ddq1 ...
    +(m2*l1*r2*cos(q2)+ m2*r2^2 + m3*l1*l2*cos(q2) + m3*l1*r3*cos(q1+q2) + m3*l1^2 + 2*m3*l2*r3*cos(q3) + m3*r3^2 + I1 + I2 + I3) .* ddq2 ...
    +(m3*l1*r3*cos(q1+q2) + m3*l2*r3*cos(q3) + m3*r3^2 + I1 + I2 + I3).*ddq3 + ((-2*m2*l1*r2*sin(q2)- 2*l1*r3*m3*sin(q2+q3) - 2*m3*l1*l2*sin(q2)).*dq1).*dq2 ...
    +((-2*m3*l1*r3*sin(q2+q3) - 2*m3*l2*r3*sin(q3)).*dq1).*dq3 + ((-m2*l1*r2*sin(q2)-m3*l1*l2*sin(q2)-m3*l1*r3*sin(q2+q3)).*dq2).*dq2 ...
    +((-2*m3*l1*r3*sin(q2+q3) - 2*m3*l2*r3*sin(q3)).*dq2).*dq3 + ((-m3*l2*r3*sin(q3)).*dq3).*dq3;

T2 = (m1*l1*r2*cos(q1)+ m2*r2^2 + m3*l1*l2*cos(q2) + m3*l1*r3*cos(q1+q2)+ m3*l1^2 + 2*m3*l2*r3*cos(q3) + m3*r3^2 + I1 + I2 + I3).*ddq1 ...
    +(m2*r2^2 + m3*l2^2 + 2*m3*l2*r3*cos(q3) +m3*r3^2 + I1 + I2 + I3).*ddq2 + (m3*l3*r3*cos(q3) + m3*r3^2 + I1 + I2 + I3).*ddq3 ...
    +((m1*l1*r2*sin(q2) + m3*l1*r3*sin(q2+q3) + m3*l1*l2*sin(q2)).*dq1).*dq1 + ((-m3*l3*r3*sin(q3)).*dq3).*dq3 + ((-2*m3*l2*r3*sin(q3)).*dq1).*dq3 ...
    +((-2*m3*l2*r3*sin(q3)).*dq2).*dq3;

T3 = (m3*l1*r3*cos(q2+q3) + l2*r3*cos(q3) + m3*r3^2 + I1 + I2 + I3).*ddq1 + (m3*l3*r3*cos(q3)+ m3*r3^2 + I1 + I2 + I3).*ddq2 + (m3*r3^2 + I1 + I2 + I3).*ddq3 ...
    + ((m3*l1*r3*sin(q2+q3) + m3*l2*r3*sin(q3)).*dq1).*dq1 + ((2*m3*l2*r3*sin(q3)).*dq1).*dq3 + ((m3*l2*r3*sin(q3)).*dq2).*dq2;

%--------------------------------------------------------------------------
%plotting joint parameter - torque
figure(1);
plot(t,T1,'b');
title('Joint 1 Torque v/s time'); xlabel('time (secs)');ylabel('T1 (N-m)');
grid on; grid minor;

figure(2);
plot(t,T2,'b');
title('Joint 2 Torque v/s time'); xlabel('time (secs)');ylabel('T2 (N-m)');
grid on; grid minor;

figure(3);
plot(t,T3,'b'); title('Joint 3 Torque v/s time');
xlabel('time (secs)');ylabel('T3 (N-m)');
grid on; grid minor;

end