function [Y_CSt,Z_CSt] = TeeConvSt(Q_s,Q_c,F_s,F_c,F_st,rho,ReC,Z_CSiL)
% TeeDivSide Computes the local loss coefficient for converter in straight
% passage
%   Z_CSt = TeeConvSt(Q_s,Q_c,F_s,F_c,F_st) = DeltaP/(rho*w_c^2*0.5)
%                   ___________________
%
%                   st ->         c ->
%                   ______       ______
%                         |  A  |
%                         |  |  |
%                         |  s  |
%   Subscript used                   Variables
%   st = straight passage            w  = fluid velocity
%   s  = side passage                Q  = volume flow rate
%   c  = combined passage            F  = cross area
%                                    rho = fluid density [km/m3]
%                                    M = mass flow rate [kg/s]
%
%   This function computes the local loss coefficient Z_CSt for a converter in
%   straight passage according to Idelchick "Handbook of Hydraulic Resistance",
%   CRC press, 3rd edition, 1994 (page 432). It also transforms the same 
%   coefficient into Y_CSt, which expresses the pressure drop as function
%   of the combined mass flow rate M_c, so that
%   Dp= Z_CSt*0.5*rho*w_c^2 = Y_CSt*M_c^2

ReL=3500;   % Re below which laminar equations are used
ReT=4000;   % Re above which turbulent equations are used
w_s=Q_s./F_s;
w_c=Q_c./F_c;

if F_s+F_st>F_c & F_st==F_c
  Z_CStT=1.55*Q_s./Q_c-(Q_s./Q_c).^2; % [-] res.coef. by Idelchick in turbulent
  Z_CStT=Z_CStT*2.2; % adding the correction introduced by Parasol report
  a0=(F_s./F_c<=0.35).*(1.8-Q_s./Q_c)+...
     (F_s./F_c>0.35).*((Q_s./Q_c<=0.2).*(1.8-4*Q_s./Q_c)+(Q_s./Q_c>0.2).*(1.2-Q_s./Q_c));
  Z_CStL=2*Z_CSiL+a0.*(1-Q_s./Q_c).^2-(1.6-0.3*F_s./F_c).*(w_s./w_c).^2;
  Z_CSt=(ReC>=ReT).*Z_CStT+...
        (ReC<=ReL).*Z_CStL+...
        (ReC>ReL).*(ReC<ReT).*(Z_CStL+(ReC-ReL)/(ReT-ReL).*(Z_CStT-Z_CStL));
  Y_CSt=Z_CSt*0.5./(F_c.^2*rho);    % [1/m.kg] Dp=Y_DSi*M_c^2
else
  warning('Cross area combination of Tee-piece not covered!')
end
end

