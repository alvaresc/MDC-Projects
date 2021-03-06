function [c,ceq] = funcNL(xVec,Mb,MG,Mc)
DB = xVec(1);
DC = xVec(2);

% ceq = rem(xVec(2),1/16);
% ceq(1) = xVec(1) - 1/16*floor(xVec(1)/(1/16));
% ceq = 1/16-DC+DB;
ceq = [];

% Sut = 148000; %psi
Sut = 85000;
% find Se
Cload = 1;
CsizeB = 0.869*DB^(-0.097);
CsizeG = CsizeB;
CsizeC = 0.869*DC^(-0.097);
% Csurf = 2.7*148^-.265;
Csurf = 2.7*85^-.265;
Ctemp = 1;
Creliab = 0.814;
Sep = 0.5*Sut;
SeB = Cload*CsizeB*Csurf*Ctemp*Creliab*Sep;
SeG = Cload*CsizeG*Csurf*Ctemp*Creliab*Sep;
SeC = Cload*CsizeC*Csurf*Ctemp*Creliab*Sep;

% find fatigue stress conectration factors at G
[Kt_bending,Kt_torsion] = interpFatigue(DC,DB);
q = 0.862;
Kf_bending = 1+q*(Kt_bending-1);
Kf_torsion = 1+q*(Kt_torsion-1);
sigma_G = Kf_bending*((32*MG)/(pi*DB^3));
% T = 600;
T = 300;
tau_torsion_G = Kf_torsion*((16*T)/(pi*DB^3));
sigma_m_G = sqrt(3)*tau_torsion_G;
sigma_m_B = 0;
sigma_m_C = (16*T*sqrt(3))/(pi*DC^3);

sigma_B = (32*Mb)/(pi*DB^3);
sigma_C = (32*Mc)/(pi*DC^3);

c(1) = -1/1.5 + sigma_m_B/Sut + sigma_B/SeB;
c(2) = -1/1.5 + sigma_m_G/Sut + sigma_G/SeG;
c(3) = -1/1.5 + sigma_m_C/Sut + sigma_C/SeC;
c(4) = DC/DB-2;
c(5) = 1.09-DC/DB;

end