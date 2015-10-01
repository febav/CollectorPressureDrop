function [Y_CSi,Z_CSi,Z_CSiL] = TeeConvSide(Q_s,Q_c,F_s,F_c,F_st,rho,ReC)
% TeeDivSide Computes the local loss coefficient for diverter in side
% passage
%   Z_CSi = TeeConvSide(Q_s,Q_c,F_s,F_c,F_st,rho) = DeltaP/(rho*w_c^2*0.5)
%                   ___________________
%
%                   st ->         c ->
%                   ______       ______
%                         |  A  |
%                         |  |  |
%                         |  s  |
%   Subscript used                   Variables
%   st = straight passage            w  = fluid velocity [m/s]
%   s  = side passage                Q  = volume flow rate [m3/s]
%   c  = combined passage            F  = cross area [m2]
%                                    M = mass flow rate [kg/s]
%
%   This function computes the local loss coefficient Z_CSi for a converter in
%   side passage according to Idelchick "Handbook of Hydraulic Resistance",
%   CRC press, 3rd edition, 1994 (page 432). It also transforms the same 
%   coefficient into Y_CSi, which expresses the pressure drop as function
%   of the combined mass flow rate M_c, so that
%   Dp= Z_CSi*0.5*rho*w_c^2 = Y_CSi*M_c^2

ReL=3500;   % Re below which laminar equations are used
ReT=4000;   % Re above which turbulent equations are used

A=zeros(size(Q_s));
w_s=Q_s./F_s;
w_c=Q_c./F_c;
wratio=w_s./w_c;
Fratio=F_s./F_c;
Qratio=Q_s./Q_c;

if F_s+F_st>F_c & F_st==F_c
    A=(Fratio<=0.35)*1+...
    (Fratio>0.35).*(Qratio<=0.4)*0.9.*(1-Qratio)+...
    (Fratio>0.35).*(Qratio>0.4)*0.55;   
else
    warning('Cross area combination of Tee-piece not covered!')
end
Z_CSiT=A.*(1+(wratio).^2-2*(1-Qratio).^2);    % [-] turbulent res.coef. by Idelchick
Z_CSiL=2*Z_CSiT+150./ReC;         % [-] laminar res.coef. by Idelchick
Z_CSi=(ReC>=ReT).*Z_CSiT+...
      (ReC<=ReL).*Z_CSiL+...
      (ReC>ReL).*(ReC<ReT).*(Z_CSiL+(ReC-ReL)/(ReT-ReL).*(Z_CSiT-Z_CSiL));

Y_CSi=Z_CSi*0.5./(F_c.^2.*rho);      % [1/m.kg] Dp=Y_CSi*M_c^2
end