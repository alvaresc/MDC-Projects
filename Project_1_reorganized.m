%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ian Landwehr, Sam Alvares, and Sam Ridgley
% March 23rd, 2021
%
% ME480: Machine Component Design
% Dr. Constans
% Spring 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear variables

%% Basic Level
% define constants
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
w = 0.712; % inches
t = 0.09 ; % since values for arm and link are the same for simplicity
Force = 20 ; % lbf (refers to half the force applied to the handle due to symmetry)
FOS = 1.8 ; % (refers to the Factor of Safety required)
alpha = 0;

% solve for reaction forces
F1 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F3 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F2_x = -(Force)*(sin(alpha)*sin(theta)*y3 - sin(alpha)*sin(theta)*y4 + sin(alpha)*cos(theta)*x3 - sin(theta)*cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % x-component of the force present at rivet 2
F2_y = -(Force)*(sin(alpha)*cos(theta)*y4 - sin(theta)*cos(alpha)*y3 - cos(alpha)*cos(theta)*x3 + cos(alpha)*cos(theta)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % y-component of the force present at rivet 2
F2 = sqrt( (F2_x).^2 + (F2_y).^2 ) ;
Flink = F1;

% print reactions
fprintf('------ Basic Level (alpha=0) ------\n')
fprintf('Reaction forces:\n')
fprintf('F1 = %4.4f\n',F1)
fprintf('F2x = %4.4f\n',F2_x)
fprintf('F2y = %4.4f\n',F2_y)
fprintf('F2 = %4.4f\n',F2)
fprintf('F3 = %4.4f\n',F3)

% rivet strength - direct shear
pin_radius_outer = 1/8 ; % inches
pin_thickness = 1/16 ; % inches
pin_radius_inner = pin_radius_outer - pin_thickness ; % inches
pin_inner_area = pi * pin_radius_inner^2 ; % inches^2
pin_outer_area = pi * pin_radius_outer^2 ; % inches^2
Area_Pin =  pin_outer_area - pin_inner_area;  % inches^2
Tau_Pin_1_Shear = F1 / Area_Pin;  % psi 
Tau_Pin_2_Shear = F2 / Area_Pin; % psi
Tau_Pin_3_Shear = F3 / Area_Pin;  % psi
tau_rivet_critical = Tau_Pin_1_Shear;
Sy_yield_rivet = FOS*2*tau_rivet_critical;

% link strength - normal stress (axial)
Aaxial = t*(w-pin_diameter);
sigma_axial = Flink/Aaxial;
Sy_yield_axial = FOS*sigma_axial;

% link strenght - normal stress (bearing)
Abearing = t*pin_diameter;
sigma_bearing = Flink/Abearing;
Sy_yield_bearing = FOS*sigma_bearing;

% link strength
l = sqrt((w/2)^2-(pin_diameter/2)^2); 
Atearout = 2*l*t;
tau_tearout = Flink/Atearout;
Sy_yield_tearout = FOS*2*tau_tearout;

% print yield strengths
fprintf('\nYield strengths:\n')
fprintf('Sy_yield_rivet = %4.4f\n',Sy_yield_rivet)
fprintf('Sy_yield_axial = %4.4f\n',Sy_yield_axial)
fprintf('Sy_yield_bearing = %4.4f\n',Sy_yield_bearing)
fprintf('Sy_yield_tearout = %4.4f\n',Sy_yield_tearout)

%% Medium Level
fprintf('\n------ Medium Level (alpha = 10 deg) ------\n')
alpha = 10*pi/180; %rad

% solve for reaction forces
F1 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F3 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F2_x = -(Force)*(sin(alpha)*sin(theta)*y3 - sin(alpha)*sin(theta)*y4 + sin(alpha)*cos(theta)*x3 - sin(theta)*cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % x-component of the force present at rivet 2
F2_y = -(Force)*(sin(alpha)*cos(theta)*y4 - sin(theta)*cos(alpha)*y3 - cos(alpha)*cos(theta)*x3 + cos(alpha)*cos(theta)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % y-component of the force present at rivet 2
F2 = sqrt( (F2_x).^2 + (F2_y).^2 ) ;
Flink = F1;

% print reactions
fprintf('Reaction forces:\n')
fprintf('F1 = %4.4f\n',F1)
fprintf('F2x = %4.4f\n',F2_x)
fprintf('F2y = %4.4f\n',F2_y)
fprintf('F2 = %4.4f\n',F2)
fprintf('F3 = %4.4f\n',F3)

% rivet strength - direct shear
pin_radius_outer = 1/8 ; % inches
pin_thickness = 1/16 ; % inches
pin_radius_inner = pin_radius_outer - pin_thickness ; % inches
pin_inner_area = pi * pin_radius_inner^2 ; % inches^2
pin_outer_area = pi * pin_radius_outer^2 ; % inches^2
Area_Pin =  pin_outer_area - pin_inner_area;  % inches^2
Tau_Pin_1_Shear = F1 / Area_Pin ; % psi 
Tau_Pin_2_Shear = F2 / Area_Pin ; % psi
Tau_Pin_3_Shear = F3 / Area_Pin ; % psi
tau_rivet_critical = Tau_Pin_1_Shear;
Sy_yield_rivet = FOS*2*tau_rivet_critical;

% link strength - normal stress (axial)
Aaxial = t*(w-pin_diameter);
sigma_axial = Flink/Aaxial;
Sy_yield_axial = FOS*sigma_axial;

% link strenght - normal stress (bearing)
Abearing = t*pin_diameter;
sigma_bearing = Flink/Abearing;
Sy_yield_bearing = FOS*sigma_bearing;

% link strength
l = sqrt((w/2)^2-(pin_diameter/2)^2); 
Atearout = 2*l*t;
tau_tearout = Flink/Atearout;
Sy_yield_tearout = FOS*2*tau_tearout;

% print yield strengths
fprintf('\nYield strengths:\n')
fprintf('Sy_yield_rivet = %4.4f\n',Sy_yield_rivet)
fprintf('Sy_yield_axial = %4.4f\n',Sy_yield_axial)
fprintf('Sy_yield_bearing = %4.4f\n',Sy_yield_bearing)
fprintf('Sy_yield_tearout = %4.4f\n',Sy_yield_tearout)

%% Advanced Level
% fprintf('\n------ Advanced Level (variable alpha) ------\n')
alpha_deg = (-15:75).';
alpha = alpha_deg*(pi/180); %rad

% solve for reaction forces
F1 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F3 = Force * (sin(alpha)*y4 + cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; 
F2_x = -(Force)*(sin(alpha)*sin(theta)*y3 - sin(alpha)*sin(theta)*y4 + sin(alpha)*cos(theta)*x3 - sin(theta)*cos(alpha)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % x-component of the force present at rivet 2
F2_y = -(Force)*(sin(alpha)*cos(theta)*y4 - sin(theta)*cos(alpha)*y3 - cos(alpha)*cos(theta)*x3 + cos(alpha)*cos(theta)*x4) / (sin(theta)*y3 + cos(theta)*x3) ; % y-component of the force present at rivet 2
F2 = sqrt( (F2_x).^2 + (F2_y).^2 ) ;
Flink = F1;

% rivet strength - direct shear
pin_radius_outer = 1/8 ; % inches
pin_thickness = 1/16 ; % inches
pin_radius_inner = pin_radius_outer - pin_thickness ; % inches
pin_inner_area = pi * pin_radius_inner^2 ; % inches^2
pin_outer_area = pi * pin_radius_outer^2 ; % inches^2
Area_Pin =  pin_outer_area - pin_inner_area;  % inches^2
Tau_Pin_1_Shear = F1 / Area_Pin ; % psi 
Tau_Pin_2_Shear = F2 / Area_Pin ; % psi
Tau_Pin_3_Shear = F3 / Area_Pin ; % psi
tau_rivet_critical = Tau_Pin_1_Shear;
Sy_yield_rivet = FOS*2*tau_rivet_critical;

% link strength - normal stress (axial)
Aaxial = t*(w-pin_diameter);
sigma_axial = Flink/Aaxial;
Sy_yield_axial = FOS*sigma_axial;

% link strenght - normal stress (bearing)
Abearing = t*pin_diameter;
sigma_bearing = Flink/Abearing;
Sy_yield_bearing = FOS*sigma_bearing;

% link strength
Rlink = w/2
Rrivet = pin_diameter/2
l = sqrt((Rlink)^2-(Rrivet)^2); 
Atearout = 2*l*t;
tau_tearout = Flink/Atearout;
Sy_yield_tearout = FOS*2*tau_tearout;

% plot required strengths 
figure(1)
plot(alpha_deg, Sy_yield_rivet, '-r') ;
hold on 
plot(alpha_deg, Sy_yield_axial, '-b') ;
hold on 
plot(alpha_deg, Sy_yield_bearing, '-g') ;
hold on 
plot(alpha_deg, Sy_yield_tearout, '-c') ;
xlabel('Alpha [degrees]') ;
ylabel('Required Yield Strength [psi]') ;
legend('S_{y,rivet}', 'S_{y,axial}', 'S_{y,bearing}','S_{y,tearout})','Location','NORTHWEST')

% plot crushing force
figure(2)
F_can = 2*F1.*cos(theta);
plot(alpha_deg,F_can,'.k');
xlabel('Alpha [degrees]')
ylabel('Force on Can [lb_{f}]')

% claculate max crushing forces angle
[maxF,I] = max(F_can);
alphaMaxF = alpha_deg(I)
