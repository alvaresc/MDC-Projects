%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sam Alvares, Ian Landwehr, Sam Ridgely
% ME480-03: Machine Component Design
% Dr. Constans
% 
% Project III: Shaft Deflection in a Gear Reducer
% Due: April 6th, 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear variables; close all;
%% Given Measurements from Project Problem Statement

% Given Measurements            % Units        
F = 150;                        % lbf        
x_1 = 1.5625;                    % in            
x_2 = 2.375;                    % in
x_3 = 3.125;                    % in
x_4 = 0.6875;                   % in
D = 6;                          % in
d = 2;                          % in
L_1 = 2;                        % in
N = 1001;                       % (-)
L_output = 3.125;               % in
L_input = 2.375;                % in
E = 30*10^6;                    % psi
d_shaft = 0.5;                  % in

% Calculating EI for output and input
EI = E*pi*d_shaft^4/64;

% Solving for reaction forces
B_2 = F*x_3/x_2;
B_1 = F-B_2;
G_1 = F*L_1/(D/2);
G_2 = -G_1*(d/2)/(D/2);    % correct
B_3 = 0;
B_4 = 0;
B_6 = -(G_1*x_1)/x_2;
B_5 = -G_1-B_6;
B_8 = (G_2*x_4-G_1*x_1)/x_2;       %
B_7 = G_2-G_1-B_8; 

% Output Shaft y direction
C_1_1 = (-B_1/6)*(x_2^2);
C_2_1 = 0;

% Output Shaft x direction
C_1_2 = (G_1*D/4)*(x_2-x_1)^2/x_2;
C_2_2 = 0;

% Output Shaft z direction
C_1_3 = ((-B_5/6)*(x_2^3)-(G_1/6)*((x_2-x_1)^3))/(-x_2);
C_2_3 = 0;

% Input Shaft y direction
C_1_4 = 0;
C_2_4 = 0;

% Input Shaft x direction
C_1_5 = ((G_2*D/4)*(x_2-x_4)^2-(G_1*d/4)*(x_2-x_1)^2)/x_2;
C_2_5 = 0;

% Input Shaft z direction
C_1_6 = ((-B_7/6)*(x_2^3)+(G_2/6)*(x_2-x_4)^3-(G_1/6)*(x_2-x_1)^3)/(-x_2);
C_2_6 = 0;

% Equations for Output shaft

% y direction
[x_1_p,v_1_p,M_1_p,theta_1_p,y_1_p] = deal(zeros(N,1));
for i = 1:N
  x_1_p(i) = (i-1)*L_output/(N-1);
  v_1_p(i) = (B_1)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^0+(B_2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^0-(F)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^0;
  M_1_p(i) = (B_1)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^1+(B_2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^1-(F)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^1;
  theta_1_p(i) = ((B_1/2)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^2+(B_2/2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^2-(F/2)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^2+C_1_1)/EI;
  y_1_p(i) = ((B_1/6)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^3+(B_2/6)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^3-(F/6)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^3+C_1_1*x_1_p(i)+C_2_1)/EI;
end

figure(1)
subplot(4,1,1), plot(x_1_p,v_1_p,'LineWidth',2), grid on
title('Output shaft y direction','FontSize',20)
ylabel('Shear Force (lbf)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_1_p,M_1_p,'LineWidth',2), grid on
ylabel('Bending Moment (lbf in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,3), plot(x_1_p,theta_1_p,'LineWidth',2), grid on
ylabel('Theta (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_1_p,y_1_p,'LineWidth',2), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)

% z direction
[x_3_p,v_3_p,M_3_p,theta_3_p,y_3_p] = deal(zeros(N,1));
for i = 1:N
  x_3_p(i) = (i-1)*L_output/(N-1);
  v_3_p(i) = -(B_5)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^0-(G_1)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^0-(B_6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^0;
  M_3_p(i) = -(B_5)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^1-(G_1)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^1-(B_6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^1;
  theta_3_p(i) = (-(B_5/2)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^2-(G_1/2)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^2-(B_6/2)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^2+C_1_3)/EI;
  y_3_p(i) = (-(B_5/6)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^3-(G_1/6)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^3-(B_6/6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^3+C_1_3*x_3_p(i)+C_2_3)/EI;
end

figure(2)
subplot(4,1,1), plot(x_3_p,v_3_p,'LineWidth',2), grid on
title('Output shaft z direction','FontSize',20)
ylabel('Shear Force (lbf)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_3_p,M_3_p,'LineWidth',2), grid on
ylabel('Bending Moment (lbf in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,3), plot(x_3_p,theta_3_p,'LineWidth',2), grid on
ylabel('Theta (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_3_p,y_3_p,'LineWidth',2), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)

% Find Max Deflections
[yMaxOutput,IndexYMaxOutput] = max(abs(y_1_p))
zMaxOutput = max(abs(y_3_p))


% Equations for Input shaft

% y direction
[x_4_p,v_4_p,M_4_p,theta_4_p,y_4_p] = deal(zeros(N,1));
for i = 1:N
  x_4_p(i) = (i-1)*L_input/(N-1);
  v_4_p(i) = (B_3)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^0+(B_4)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^0;
  M_4_p(i) = (B_3)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^1+(B_4)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^1;
  theta_4_p(i) = ((B_3/2)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^2+(B_4/2)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^2+C_1_4)/EI;
  y_4_p(i) = ((B_3/6)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^3+(B_4/6)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^3+C_1_4*x_4_p(i)+C_2_4)/EI;
end

figure(3)
subplot(4,1,1), plot(x_4_p,v_4_p,'LineWidth',2), grid on
title('Input shaft y direction','FontSize',20)
ylabel('Shear Force (lbf)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_4_p,M_4_p,'LineWidth',2), grid on
ylabel('Bending Moment (lbf in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,3), plot(x_4_p,theta_4_p,'LineWidth',2), grid on
ylabel('Theta (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_4_p,y_4_p,'LineWidth',2), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)

% z direction
[x_6_p,v_6_p,M_6_p,theta_6_p,y_6_p] = deal(zeros(N,1));
for i = 1:N
  x_6_p(i) = (i-1)*L_input/(N-1);
  v_6_p(i) = -(B_7)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^0+(G_2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^0-(G_1)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^0-(B_8)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^0;
  M_6_p(i) = -(B_7)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^1+(G_2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^1-(G_1)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^1-(B_8)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^1;
  theta_6_p(i) = (-(B_7/2)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^2+(G_2/2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^2-(G_1/2)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^2-(B_8/2)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^2+C_1_6)/EI;
  y_6_p(i) = (-(B_7/6)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^3+(G_2/6)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^3-(G_1/6)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^3-(B_8/6)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^3+C_1_6*x_6_p(i)+C_2_6)/EI;
end

figure(4)
subplot(4,1,1), plot(x_6_p,v_6_p,'LineWidth',2), grid on
title('Input shaft z direction','FontSize',20)
ylabel('Shear Force (lbf)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_6_p,M_6_p,'LineWidth',2), grid on
ylabel('Bending Moment (lbf in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,3), plot(x_6_p,theta_6_p,'LineWidth',2), grid on
ylabel('Theta (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_6_p,y_6_p,'LineWidth',2), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)


% Find Max Deflections
yMaxInput = max(abs(y_1_p))
zMaxInput = max(abs(y_3_p))