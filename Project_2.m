%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ian Landwehr, Sam Alvares, and Sam Ridgley
% March 30rd, 2021
%
% ME480: Machine Component Design
% Dr. Constans
% Spring 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic Level 
clc
close all
clear variables

% Geometry
Aloc = [1 -24 36];
COMloc = [38 -36 36];
MomentArm = Aloc-COMloc;

% Find External Moments/Forces
ExtForces = [0 -30 0]; % Shear Force on A
ExtMoment = cross(MomentArm, ExtForces);

% Beam Shear due to Fy
T = 1/4 ; %     The cross sectional wall thickness area. 
Do = 1 ;%       The outer Diameter of the lower pole
Di = 1-1/4;%    The inner Diameter of the lower pole
alpha = 90; %   Degrees (half-circle)
Q = 2*sind(alpha)/3*((Do/2)^3-(Di/2)^3); 
I = pi/64*(Do^4-Di^4);
V = abs(ExtForces(2));
Tau = V * Q /I/T; 

% Bending Moment due to Mz
Mz = abs(ExtMoment(3));
sigma = ExtMoment * Do / 2 / I ;

%% Medium Level
clc
close all
clear variables

% Geometry
Bloc = [0 -180 1];
COMloc = [38 -36 36];
MomentArm = Bloc-COMloc;

% Find External Moments/Forces
ExtForces = [0 -30 0]; % Shear Force on A
ExtMoment = cross(MomentArm, ExtForces);

% Beam Shear due to Fy
Do = 2 ;%       The outer Diameter of the lower pole
Di = 1.5;%    The inner Diameter of the lower pole
alpha = 90; %   Degrees (half-circle)
Q = 2*sind(alpha)/3*((Do/2)^3-(Di/2)^3); 
I = pi/64*(Do^4-Di^4);
T = 0.5;
V = abs(ExtForces(2));
TauBeam = V * Q /I/T;

% Torsion due to Mz
T = abs(ExtMoment(3));
J = pi/32*(Do^4-Di^4);
TauTorsion = (T*Do)/(2*J)

% Bending Moment due to Mx
Mx = abs(ExtMoment(1));
sigmaMx = Mx * Do / 2 / I 

% Enter stress matricies
sigma1 = [         0, 0,TauTorsion;
                   0, 0,         0;
          TauTorsion, 0,   sigmaMx];
      
sigma2 = [0,                     0,                     0;
          0,                     0, -(TauTorsion-TauBeam);
          0, -(TauTorsion-TauBeam),                     0];

sigma3 = [         0, 0,-TauTorsion;
                   0, 0,          0;
          -TauTorsion, 0,  -sigmaMx];
      
sigma4 = [0,                     0,                    0;
          0,                     0, (TauTorsion-TauBeam);
          0,  (TauTorsion-TauBeam),                    0];
      
% Find Principle stresses
[V1,D1] = eigs(sigma1)
[V2,D2] = eigs(sigma2)
[V3,D3] = eigs(sigma3)
[V4,D4] = eigs(sigma4)











