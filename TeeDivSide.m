function [Y_DSi,Z_DSi] = TeeDivSide(Q_s,Q_c,F_s,F_c,F_st,rho,ReC)
% TeeDivSide Computes the local loss coefficient for diverter in side
% passage
%   Z_DSi = TeeDivSide(Q_s,Q_c,F_s,F_c,F_st,rho) = DeltaP/(rho*w_c^2*0.5)
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
%   This function computes the local loss coefficient Z_DSi for a diverter
%   in side passage according to Idelchick "Handbook of Hydraulic Resistance",
%   CRC press, 3rd edition, 1994 (page 418). It also transforms the same 
%   coefficient into Y_DSi, which expresses the pressure drop as function
%   of the combined mass flow rate M_c, so that
%   Dp= Z_DSi*0.5*rho*w_c^2 = Y_DSi*M_c^2

ReL=3500;   % Re below which laminar equations are used
ReT=4000;   % Re above which turbulent equations are used

Aprime=zeros(size(Q_s));
Kprime_s=Aprime;
w_s=Q_s./F_s;
w_c=Q_c./F_c;
wratio=w_s./w_c;
Qratio=Q_s./Q_c;
Fratio=F_s./F_c;

if F_s+F_st>F_c & F_st==F_c     
        Aprime=(Fratio<=0.35).*(Qratio<=0.4).*(1.1-0.7*Qratio)+...
            (Fratio<=0.35).*(Qratio>0.4)*0.85+...
            (Fratio>0.35).*(Qratio<=0.6).*(1-0.6*Qratio)+...
            (Fratio>0.35).*(Qratio>0.6)*0.6;
        Kprime_s=0;
        %     elseif F_s+F_st==F_c
        %         Aprime(i)=1;
        %         Kprime_s(i)=1;
else
    warning('Cross area combination of Tee-piece not covered!')
end
Z_DSiT=Aprime.*(1+wratio.^2)-Kprime_s.*wratio.^2; % [-] turbulent res.coef. by Idelchick
k1=(Qratio<=0.6).*(0.9+Qratio)+(Qratio>0.6).*(1.5-(Qratio-0.6)/2);
Z_DSiL=(k1+1).*Z_DSiT+150./ReC; % [-] laminar res.coef. by Idelchick
Z_DSiL=Z_DSiL*1.11; % correction from Parasol report
Z_DSi=(ReC>=ReT).*Z_DSiT+...
      (ReC<=ReL).*Z_DSiL+...
      (ReC>ReL).*(ReC<ReT).*(Z_DSiL+(ReC-ReL)/(ReT-ReL).*(Z_DSiT-Z_DSiL));
  
Y_DSi=Z_DSi*0.5./(F_c.^2.*rho);             % [1/m.kg] Dp=Y_DSi*M_c^2
end

