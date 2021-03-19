%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ian Landwehr, Sam Alvares, and Sam Ridgley
% March 23rd, 2021
%
% ME480: Machine Component Design
% Dr. Constans
% Spring 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc ; clear all ;close all ;





%% Given Measurements for the Free Body Diagram
x1 = 1.032 ;  % inches
x2 = 0 ;      % inches
x3 = 2.825 ;  % inches
x4 = 9.625 ;  % inches
y1 = 5.936 ;  % inches
y2 = 0 ;      % inches
y3 = 1.16 ;   % inches
y4 = 7.9375 ; % inches
pin_diameter = 0.25 ; % inches
theta = 20.577 * (pi/180) ; % (refers to the theta (in radians) calculated using Pythagorean theorem)
t_arm = 0.09 ; % inches (refers to the thickness of the arm)
t_link = 0.09 ; % inches (refers to the thickness of the link)
thickness = 0.09 ; % since values for arm and link are the same for simplicity
Force = 20 ; % lbf (refers to half the force applied to the handle due to symmetry)
FOS = 1.8 ; % (refers to the Factor of Safety required)
alpha_degrees = -15:75 ; % (refers to alpha angle range in degrees)
alpha = alpha_degrees * (pi/180) ; % (refers to alpha angle converted from degrees to radians)
zero_degree_index = 16 ; % (refers to the index at which alpha is zero degrees)
ten_degree_index = 26  ; % (refers to the index at which alpha is ten degrees)

%% Determining Numerical Values for Reaction Forces at Pins 1, 2, and 3

% The intermediate, variable solutions were documented in Maple using the
% system of four equations to solve for the four unknowns. The four
% reaction force values were documented as follows:

F1 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F3 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F2_x = -(Force)*(sin(alpha)*sin(theta)*y3 - sin(alpha)*sin(theta)*y4 + sin(alpha)*cos(theta)*x3 - sin(theta)*cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % x-component of the force present at rivet 2
F2_y = -(Force)*(sin(alpha)*cos(theta)*y4 - sin(theta)*cos(alpha)*y3 - cos(alpha)*cos(theta)*x3 + cos(alpha)*cos(theta)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % y-component of the force present at rivet 2
F2 = sqrt( (F2_x).^2 + (F2_y).^2 ) ; % Combine the x and y component force vectors for rivet 2. 

% Below are intermediate calculations

fprintf('This is F1 with alpha at zero degrees')
F1(zero_degree_index)
fprintf('This is F1 with alpha at ten degrees')
F1(ten_degree_index)
fprintf('This is F2x with alpha at zero degrees')
F2_x(zero_degree_index)
fprintf('This is F2x with alpha at ten degrees')
F2_x(ten_degree_index)
fprintf('This is F2y with alpha at zero degrees')
F2_y(zero_degree_index)
fprintf('This is F2y with alpha at ten degrees')
F2_y(ten_degree_index)
fprintf('This is F2 with alpha at zero degrees')
F2(zero_degree_index)
fprintf('This is F2 with alpha at ten degrees')
F2(ten_degree_index)
fprintf('This is F3 with alpha at zero degrees')
F3(zero_degree_index)
fprintf('This is F3 with alpha at ten degrees')
F3(ten_degree_index)

%% Calculation of Shear Stress (DOUBLE-CHECK THIS SECTION)

pin_radius_outer = 1/8 ; % inches
pin_thickness = 1/16 ; % inches
pin_radius_inner = pin_radius_outer - pin_thickness ; % inches
pin_inner_area = pi * pin_radius_inner^2 ; % inches^2
pin_outer_area = pi * pin_radius_outer^2 ; % inches^2


Area_Pin =  pin_outer_area - pin_inner_area ; % inches^2
Tau_Pin_1_Shear = F1 / Area_Pin ; % psi
Tau_Pin_2_Shear = F2 / Area_Pin ; % psi
Tau_Pin_3_Shear = F3 / Area_Pin ; % psi

%% Plotting Different Shear Stresses at Rivets 1, 2, and 3

figure(1)
plot(alpha_degrees, Tau_Pin_1_Shear, '-r') ;
hold on 
plot(alpha_degrees, Tau_Pin_2_Shear, '-b') ;
hold on 
plot(alpha_degrees, Tau_Pin_3_Shear, 'ko') ;
xlabel('Alpha (degrees)') ;
ylabel('Shear Stress (psi)') ;
title(' Tau (Shear) versus Angle for Rivets 1, 2, and 3') ;
legend('Rivet 1', 'Rivet 2', 'Rivet 3') ;



%% Maximum Shear Strength Calculation

% Making an array with the maximum shear stress of any pin occurs at Rivets
% 1/3 as shown in Figure 1. We therefore know that the maximum shear stress
% is always on pins located at rivet 1 and rivet 3 (with the same value). 

Maximum_Tau_Shear = Tau_Pin_3_Shear ;

% Calculate the shear strength given the FOS above and the assumption that
% a material's strength in shear is 1/2 its strength in tension. 

Shear_Strength = 2 * Maximum_Tau_Shear * FOS ;

%% Calculation of Bearing Stress

% We know that the bearing stress will be the same in the link and the arm
% beacuse they have the same diameter and thickness. 

Area_Bearing = pin_diameter * thickness ;
Sigma_Pin_1_Bearing = F1 / Area_Bearing ;
Sigma_Pin_2_Bearing = F2 / Area_Bearing ;
Sigma_Pin_3_Bearing = F3 / Area_Bearing ;

%% Plotting Different Bearing Stresses at Rivets 1, 2, and 3
figure(2)
plot(alpha_degrees, Sigma_Pin_1_Bearing, '-r') ;
hold on 
plot(alpha_degrees, Sigma_Pin_2_Bearing, '-b') ;
hold on 
plot(alpha_degrees, Sigma_Pin_3_Bearing, 'ko') ;
xlabel('Alpha (degrees)') ;
ylabel('Sigma (psi)') ;
title('Sigma (Bearing) versus Angle for Rivets 1, 2, and 3') ;
legend('Rivet 1', 'Rivet 2', 'Rivet 3') ;



%% Maximum Bearing Strength Calculation

% From Figure 2, we know the highest bearing stress will occur at pins 1
% and 3 because they experience the highest force and all the pins have the
% same area

Maximum_Sigma_Bearing = Sigma_Pin_3_Bearing ;

% Calculate the bearing strength given the FOS above

Bearing_Strength = Maximum_Sigma_Bearing * FOS ;



%% Calculation of the Tearout Stress

x = sqrt((.42)^2-(.125)^2) ;
Area_Tearout = x * thickness ;

Tau_Pin_1_Tearout = F1 / Area_Tearout ;
Tau_Pin_2_Tearout = F2 / Area_Tearout ;
Tau_Pin_3_Tearout = F3 / Area_Tearout ;



%% Plotting Different Tearout Stresses at Rivets 1, 2, and 3

figure(3)
plot(alpha_degrees, Tau_Pin_1_Tearout, '-r') ;
hold on 
plot(alpha_degrees, Tau_Pin_2_Tearout, '-b') ;
hold on 
plot(alpha_degrees, Tau_Pin_3_Tearout, 'ko') ;
xlabel('Alpha (degrees)') ;
ylabel('Tau (psi)') ;
title('Tau (Tearout) versus Angle for Rivets 1, 2, and 3') ;
legend('Rivet 1', 'Rivet 2', 'Rivet 3') ;



%% Maximum Tearout Strength Calculation 

% From Figure 3, we know the highest tearout stress will occur at pins 1
% and 3 because they experience the highest force and all the pins have the
% same area. 

Maximum_Tau_Tearout = Tau_Pin_3_Tearout ;

% Calculate the tearout strength given the FOS above and the assumption
% that a material's strength in shear is 1/2 its strength in tension.

Tearout_Strength = 2* Maximum_Tau_Tearout * FOS ;


%% Calculation of the Axial Stress

Area_Axial = thickness * .462 ;

Sigma_Pin_1_Axial = F1 / Area_Axial ;
Sigma_Pin_2_Axial = F2 / Area_Axial ;
Sigma_Pin_3_Axial = F3 / Area_Axial ;



%% Plotting Different Axial Stresses in the Link

figure(4)
plot(alpha_degrees, Sigma_Pin_1_Axial, '-r') ;
hold on 
plot(alpha_degrees, Sigma_Pin_2_Axial, '-b') ;
hold on 
plot(alpha_degrees, Sigma_Pin_3_Axial, 'ko') ;
xlabel('Alpha (degrees)') ;
ylabel('Sigma (psi)') ;
title('Sigma (Axial) versus Angle in the Link') ;
legend('Rivet 1', 'Rivet 2', 'Rivet 3') ;



%% Maximum Axial Strength Calculation

Maximum_Sigma_Axial = Sigma_Pin_3_Axial ;

% Calculate the axial strength given the FOS above

Axial_Strength = Maximum_Sigma_Axial * FOS ;



%% Completing the "Basic" Level Calculations (alpha is 0 degrees)

Shear_Strength_0 = Shear_Strength(zero_degree_index) ;
Bearing_Strength_0 = Bearing_Strength(zero_degree_index) ;
Tearout_Strength_0 = Tearout_Strength(zero_degree_index) ;
Axial_Strength_0 = Axial_Strength(zero_degree_index) ;



%% Completing the "Medium" Level Calculations (alpha is 10 degrees)

Shear_Strength_10 = Shear_Strength(ten_degree_index) ;
Bearing_Strength_10 = Bearing_Strength(ten_degree_index) ;
Tearout_Strength_10 = Tearout_Strength(ten_degree_index) ;
Axial_Strength_10 = Axial_Strength(ten_degree_index) ;



%% Completing the "Advanced" Level Calculations
% Computing Force on Can

% From the FBD the force on the can can be calculated from the force on the
% pin C
figure(5)
F_can = 2*F1.*cos(theta);
plot(alpha_degrees,F_can,'.k');
title('Force on the Can versus Angle')
xlabel('Angle (degrees)')
ylabel('Force (lbf)')


figure(6)
plot(alpha_degrees, Shear_Strength,'.m') ;
hold on
plot(alpha_degrees, Bearing_Strength,'.b') ;
hold on
plot(alpha_degrees, Tearout_Strength,'.g') ;
hold on
plot(alpha_degrees, Axial_Strength, '.r') ;
hold off

title('Strengths with Respect to Stress Type') ;
xlabel('Alpha (Degrees)') ;
ylabel('Stress (psi)') ;
legend('Shear', 'Bearing', 'Tearout', 'Axial') ;