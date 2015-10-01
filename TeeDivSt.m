function [Y_DSt,Z_DSt] = TeeDivSt(Q_s,Q_c,F_s,F_c,F_st,rho,ReC)
% TeeDivSide Computes the local loss coefficient for diverter in straight
% passage
%   Z_DSt = TeeDivSt(Q_s,Q_c,F_s,F_c,F_st,rho) = DeltaP/(rho*w_c^2*0.5)
%                   ___________________
%
%                   c ->         st ->
%                   ______       ______
%                         |  |  |
%                         |  V  |
%                         |  s  |
%   Subscript used                   Variables
%   st = straight passage            w  = fluid velocity [m/s]
%   s  = side passage                Q  = volume flow rate [m3/s]
%   c  = combined passage            F  = cross area [m2]
%                                    rho = fluid density [km/m3]
%                                    M = mass flow rate [kg/s]
%
%   This function computes the local loss coefficient Z_DSt for a diverter in
%   straight passage according to Idelchick "Handbook of Hydraulic Resistance",
%   CRC press, 3rd edition, 1994 (page 419+453 diagram 7-20). It also transforms
%   the same coefficient into Y_CSt, which expresses the pressure drop as function
%   of the combined mass flow rate M_c, so that
%   Dp= Z_DSi*0.5*rho*w_c^2 = Y_DSi*M_c^2

% w_s=Q_s/F_s;
% w_c=Q_c/F_c;

ReL=3500;   % Re below which laminar equations are used
ReT=4000;   % Re above which turbulent equations are used

tau_st=zeros(size(Q_s));
Fratio=F_s./F_c;
Qratio=Q_s./Q_c;

if F_s+F_st>F_c & F_st==F_c
  tau_st=(Fratio<=0.4)*0.4+...
         (Fratio>0.4).*(Qratio<=0.5)*2.*(2*Qratio-1)+...
         (Fratio>0.4).*(Qratio>0.5)*0.3.*(2*Qratio-1);
else
  warning('Cross area combination of Tee-piece not covered!')
end
Z_DStT=tau_st.*Qratio.^2;     % [-] turbulent res.coef. by Idelchick
Z_DStL=3*Z_DStT+33./ReC; % [-] laminar res.coef. by Idelchick
Z_DSt=(ReC>=ReT).*Z_DStT+...
      (ReC<=ReL).*Z_DStL+...
      (ReC>ReL).*(ReC<ReT).*(Z_DStL+(ReC-ReL)/(ReT-ReL).*(Z_DStT-Z_DStL));
Y_DSt=Z_DSt*0.5./(F_c.^2.*rho);  % [1/m.kg] Dp=Y_Dst*M_c^2
end

