function [Kt_bending,Kt_torsion] = interpFatigue(D,d)

DdTorsionVec = [2,1.33,1.2,1.09];
ATorsionVec = [0.86331, 0.84897, 0.83425, 0.90337];
bTorsionVec = [-.23865, -.23161, -.21649, -.12692];
DdBendingVec = [6,3,2,1.5,1.2,1.1,1.07,1.05,1.03,1.02,1.01];
ABendingVec = [0.87868,0.89334,0.90879,0.93836,0.97098,0.95120,0.97527,0.98137,0.98061,0.96048,0.91938];
bBendingVec = [-.33243,-.30860,-.28598,-.25759,-.21796,-.23757,-.20958,-.19653,-.18381,-.17711,-.17032];
r = 0.05;
rd = r/d;
Dd = [D/d];

ATorsion = interp1(DdTorsionVec,ATorsionVec,Dd);
bTorsion = interp1(DdTorsionVec,bTorsionVec,Dd);
ABending = interp1(DdBendingVec,ABendingVec,Dd);
bBending = interp1(DdBendingVec,bBendingVec,Dd);

Kt_torsion = ATorsion*rd^bTorsion;
Kt_bending = ABending*rd^bBending;

% if isnan(Kt_torsion)
%     Kt_torsion = Kt_torsion*;
%     
% elseif isnan(Kt_bending)
%     Kt_bending = 10000;
end
    
