%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sam Alvares, Ian Landwehr, Sam Ridgely
% ME480-03: Machine Component Design
% Dr. Constans
% 
% Project III: Shaft Deflection in a Gear Reducer
% Due: April 6th, 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear variables; close all;
%% Given Measurements from Project Problem Statement

% Given Measurements            % Units        
F = 150;                        % lbf        
x_1 = 1.5625;                   % in            
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

% Enter reaction forces (From Maple)
B_1 = -900/19;
B_2 = 3750/19; 
B_3 = 0;
B_4 = 0; 
B_5 = -650/19;
B_6 = -1250/19;
B_7 = 1100/19;
B_8 = 4300/57;
G_1 = 100; 
G_2 = 100/3;

% Output Shaft y direction
C_1_1 = (-B_1/6)*(x_2^2);
C_2_1 = 0;

% Output Shaft z direction
C_1_3 = ((B_5/6)*(x_2^3)+(G_1/6)*((x_2-x_1)^3))/(-x_2);
C_2_3 = 0;

% Input Shaft y direction
C_1_4 = 0;
C_2_4 = 0;

% Input Shaft z direction
C_1_6 = ((B_7/6)*(x_2^3)-(G_2/6)*(x_2-x_4)^3-(G_1/6)*(x_2-x_1)^3)/(-x_2);
C_2_6 = 0;

%% Output shaft

% y direction
[x_1_p,v_1_p,M_1_p,theta_1_p,y_1_p] = deal(zeros(N,1));
for i = 1:N
  x_1_p(i) = (i-1)*L_output/(N-1);
  v_1_p(i) = (B_1)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^0+(B_2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^0-(F)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^0;
  M_1_p(i) = (B_1)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^1+(B_2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^1-(F)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^1;
  theta_1_p(i) = ((B_1/2)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^2+(B_2/2)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^2-(F/2)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^2+C_1_1)/EI;
  y_1_p(i) = ((B_1/6)*heaviside(x_1_p(i)-0)*(x_1_p(i)-0)^3+(B_2/6)*heaviside(x_1_p(i)-x_2)*(x_1_p(i)-x_2)^3-(F/6)*heaviside(x_1_p(i)-x_3)*(x_1_p(i)-x_3)^3+C_1_1*x_1_p(i)+C_2_1)/EI;
end

% z direction
[x_3_p,v_3_p,M_3_p,theta_3_p,y_3_p] = deal(zeros(N,1));
for i = 1:N
  x_3_p(i) = (i-1)*L_output/(N-1);
  v_3_p(i) = (B_5)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^0+(G_1)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^0+(B_6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^0;
  M_3_p(i) = (B_5)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^1+(G_1)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^1+(B_6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^1;
  theta_3_p(i) = ((B_5/2)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^2+(G_1/2)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^2+(B_6/2)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^2+C_1_3)/EI;
  y_3_p(i) = ((B_5/6)*heaviside(x_3_p(i)-0)*(x_3_p(i)-0)^3+(G_1/6)*heaviside(x_3_p(i)-x_1)*(x_3_p(i)-x_1)^3+(B_6/6)*heaviside(x_3_p(i)-x_2)*(x_3_p(i)-x_2)^3+C_1_3*x_3_p(i)+C_2_3)/EI;
end

% Find Max Deflections and locations
[yMaxOutput,IndexYMaxOutput] = max(abs(y_1_p));
[zMaxOutput,IndexZMaxOutput] = max(abs(y_3_p));
yMaxLocationOutput = x_1_p(IndexYMaxOutput);
zMaxLocationOutput = x_3_p(IndexZMaxOutput);

% Find Max Bending Moment and location
[yMaxMomentOutput,IndexYMaxMomentOutput] = max(abs(M_1_p));
yMaxLocationMomentOutput = x_1_p(IndexYMaxMomentOutput);
[zMaxMomentOutput,IndexZMaxMomentOutput] = max(abs(M_3_p));
zMaxLocationMomentOutput = x_3_p(IndexZMaxMomentOutput);

% Plot

figure(1)
subplot(4,1,1), plot(x_1_p,v_1_p,'g','LineWidth',1), grid on
title('Output Shaft Y-Direction','FontSize',20)
ylabel('Shear Force (lb_{f})','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_1_p,M_1_p,'g','LineWidth',1), grid on
ylabel('Bending Moment (lb_{f}-in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(yMaxLocationMomentOutput,M_1_p(IndexYMaxMomentOutput),'ro')
label = sprintf('(%1.3f,%3.1f)',yMaxLocationMomentOutput,M_1_p(IndexYMaxMomentOutput));
text(yMaxLocationMomentOutput,M_1_p(IndexYMaxMomentOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
subplot(4,1,3), plot(x_1_p,theta_1_p,'g','LineWidth',1), grid on
ylabel('Slope (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_1_p,y_1_p,'g','LineWidth',1), grid on
hold on
plot(yMaxLocationOutput,y_1_p(IndexYMaxOutput),'ro')
label = sprintf('(%1.3f,%1.3d)',yMaxLocationOutput,y_1_p(IndexYMaxOutput));
text(yMaxLocationOutput,y_1_p(IndexYMaxOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left');
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)

figure(2)
subplot(4,1,1), plot(x_3_p,v_3_p,'g','LineWidth',1), grid on
title('Output Shaft Z-Direction','FontSize',20)
ylabel('Shear Force (lb_{f})','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_3_p,M_3_p,'g','LineWidth',1), grid on
ylabel('Bending Moment (lb_{f}-in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(zMaxLocationMomentOutput,M_3_p(IndexZMaxMomentOutput),'ro')
label = sprintf('(%1.3f,%3.1f)',zMaxLocationMomentOutput,M_3_p(IndexZMaxMomentOutput));
text(zMaxLocationMomentOutput,M_3_p(IndexZMaxMomentOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
subplot(4,1,3), plot(x_3_p,theta_3_p,'g','LineWidth',1), grid on
ylabel('Slope (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_3_p,y_3_p,'g','LineWidth',1), grid on
hold on
plot(zMaxLocationOutput,y_3_p(IndexZMaxOutput),'ro')
label = sprintf('(%1.3f,%1.3d)',zMaxLocationOutput,y_3_p(IndexZMaxOutput));
text(zMaxLocationOutput,y_3_p(IndexZMaxOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','right');
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)

%% Input shaft

% y direction
[x_4_p,v_4_p,M_4_p,theta_4_p,y_4_p] = deal(zeros(N,1));
for i = 1:N
  x_4_p(i) = (i-1)*L_input/(N-1);
  v_4_p(i) = (B_3)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^0+(B_4)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^0;
  M_4_p(i) = (B_3)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^1+(B_4)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^1;
  theta_4_p(i) = ((B_3/2)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^2+(B_4/2)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^2+C_1_4)/EI;
  y_4_p(i) = ((B_3/6)*heaviside(x_4_p(i)-0)*(x_4_p(i)-0)^3+(B_4/6)*heaviside(x_4_p(i)-x_2)*(x_4_p(i)-x_2)^3+C_1_4*x_4_p(i)+C_2_4)/EI;
end

% z direction
[x_6_p,v_6_p,M_6_p,theta_6_p,y_6_p] = deal(zeros(N,1));
for i = 1:N
  x_6_p(i) = (i-1)*L_input/(N-1);
  v_6_p(i) = (B_7)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^0-(G_2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^0-(G_1)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^0+(B_8)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^0;
  M_6_p(i) = (B_7)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^1-(G_2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^1-(G_1)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^1+(B_8)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^1;
  theta_6_p(i) = ((B_7/2)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^2-(G_2/2)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^2-(G_1/2)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^2+(B_8/2)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^2+C_1_6)/EI;
  y_6_p(i) = ((B_7/6)*heaviside(x_6_p(i)-0)*(x_6_p(i)-0)^3-(G_2/6)*heaviside(x_6_p(i)-x_4)*(x_6_p(i)-x_4)^3-(G_1/6)*heaviside(x_6_p(i)-x_1)*(x_6_p(i)-x_1)^3+(B_8/6)*heaviside(x_6_p(i)-x_2)*(x_6_p(i)-x_2)^3+C_1_6*x_6_p(i)+C_2_6)/EI;
end

% Find Max Deflections and locations
[yMaxOutput,IndexYMaxOutput] = max(abs(y_4_p));
[zMaxOutput,IndexZMaxOutput] = max(abs(y_6_p));
yMaxLocationOutput = x_4_p(IndexYMaxOutput);
zMaxLocationOutput = x_6_p(IndexZMaxOutput);

% Find Max Bending Moment and location
[yMaxMomentOutput,IndexYMaxMomentOutput] = max(abs(M_4_p));
yMaxLocationMomentOutput = x_4_p(IndexYMaxMomentOutput);
[zMaxMomentOutput,IndexZMaxMomentOutput] = max(abs(M_6_p));
zMaxLocationMomentOutput = x_6_p(IndexZMaxMomentOutput);

% Plot
figure(3)
subplot(4,1,1), plot(x_4_p,v_4_p,'g','LineWidth',1), grid on
title('Input Shaft Y-Direction','FontSize',20)
ylabel('Shear Force (lb_{f})','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_4_p,M_4_p,'g','LineWidth',1), grid on
ylabel('Bending Moment (lb_{f}-in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(yMaxLocationMomentOutput,M_4_p(IndexYMaxMomentOutput),'ro')
label = sprintf('(%1.3f,%3.1f)',yMaxLocationMomentOutput,M_4_p(IndexYMaxMomentOutput));
text(yMaxLocationMomentOutput,M_4_p(IndexYMaxMomentOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
subplot(4,1,3), plot(x_4_p,theta_4_p,'g','LineWidth',1), grid on
ylabel('Slope (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_4_p,y_4_p,'g','LineWidth',1), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(yMaxLocationOutput,y_4_p(IndexYMaxOutput),'ro')
label = sprintf('(%1.3f,%1.3d)',yMaxLocationOutput,y_4_p(IndexYMaxOutput));
text(yMaxLocationOutput,y_4_p(IndexYMaxOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left');

figure(4)
subplot(4,1,1), plot(x_6_p,v_6_p,'g','LineWidth',1), grid on
title('Input Shaft Z-Direction','FontSize',20)
ylabel('Shear Force (lb_{f})','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,2), plot(x_6_p,M_6_p,'g','LineWidth',1), grid on
ylabel('Bending Moment (lb_{f}-in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(zMaxLocationMomentOutput,M_6_p(IndexZMaxMomentOutput),'ro')
label = sprintf('(%1.3f,%3.1f)',zMaxLocationMomentOutput,M_6_p(IndexZMaxMomentOutput));
text(zMaxLocationMomentOutput,M_6_p(IndexZMaxMomentOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
subplot(4,1,3), plot(x_6_p,theta_6_p,'g','LineWidth',1), grid on
ylabel('Slope (-)','FontSize',12)
xlabel('Position (in)','FontSize',12)
subplot(4,1,4), plot(x_6_p,y_6_p,'g','LineWidth',1), grid on
ylabel('Displacement (in)','FontSize',12)
xlabel('Position (in)','FontSize',12)
hold on
plot(zMaxLocationOutput,y_6_p(IndexZMaxOutput),'ro')
label = sprintf('(%1.3f,%1.3d)',zMaxLocationOutput,y_6_p(IndexZMaxOutput));
text(zMaxLocationOutput,y_6_p(IndexZMaxOutput),label,'VerticalAlignment','bottom','HorizontalAlignment','center');

