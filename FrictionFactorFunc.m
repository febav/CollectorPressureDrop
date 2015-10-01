function ffrV=FrictionFactorFunc_scalar(ReV,D,ep,option)

% function ffr=FrictionFactorFunc(Re,D,ep,option)
% where input parameters are:
% * ReV [-] Reynolds number (can be a vector)
% * D [m] inner diameter of the pipe
% * ep [m] roughness of the pipe
% * option [-] is a parameter deciding which set of equations is used to
%   compute the friction factor.
%   option = 1-3: piece-wise function with use of different correlations in
%   the turbulent region, 64/Re in the laminar region and linear
%   interpolation for the transition region (limit values of Re defining
%   end of laminar and beginning of turbulent need to be defined inside the
%   function script. More specifically:
%   option = 1: Blasius correlation (ffr = 0.3164/Re^0.25)
%   option = 2: Colebrook correlation (1/sqrt(ffr)=(-2*log10(ep/(3.7*D)+...
%   option = 3: Haaland correlation (ffr=(-1.8*log10((ep/(3.7*D))^1.11+...
%   option = 4: overall Joseph and Yang correlation

% limit for laminar and transition region for option = 1-3
ReL=2300;   % Re below which flow is considered laminar
ReT=3100;   % Re above which flow is considered turbulent

ffrV=zeros(1,length(ReV));

switch option
    case 1      % Blasius correlation
        PoO=1;
        HoC=3;
        if ep~=0
            error('Blasius correlation can''t be used for rough pipes.')
        end
    case 2      % Colebrook correlation
        PoO=1;
        HoC=2;
    case 3      % Haaland correlation
        PoO=1;
        HoC=1;
    case 4      % Joseph and Yang correlation
        PoO=2;
        HoC=3;
        if ep~=0
            error('Joseph&Yang eq. can''t be used for rough pipes.')
        end
    otherwise
        error('Invalid "option" input passed to function. Use the help');
end
 
for i=1:length(ReV)
Re=ReV(i);      % Re needs to be a scalar
if PoO==1
    if Re<=ReT       % Laminar (+transition)
        ffrL=64/min(Re,ReL);
        ffr=ffrL;
    end

    if Re>=ReL      % Turbulent  (+transition)
        if option>1    % rough pipes
            if HoC==1           % Haaland=1
                ffrT=(-1.8*log10((ep/(3.7*D))^1.11+6.9/max(Re,ReT)))^-2;
                ffr=ffrT;
            elseif HoC==2       % Colebrook=1
                fun_Colebrook=@(ffr) 1/sqrt(ffr)-...
                    +(-2*log10(ep/(3.7*D)+2.51/max(Re,ReT)/sqrt(ffr)));
                Colebrook_1g=(-1.8*log10((ep/(3.7*D))^1.11+...
                    6.9/max(Re,ReT)))^-2; % using Halland eq
                ffrT=fzero(fun_Colebrook,Colebrook_1g);
                ffr=ffrT;
            else
                error('Error in HoC');
            end
        elseif option==1        % Blasius (smooth pipes)
            ffrT=0.3164/max(Re,ReT)^0.25;
            ffr=ffrT;
        end
    end

    if Re>=ReL && Re<=ReT     % transition region
        ffr=ffrL+(Re-ReL)/(ReT-ReL)*(ffrT-ffrL); % linear interpolation L-T
    end
elseif PoO==2      % Joseph and Yang correlation for SMOOTH pipes
    Rec=2900;
    Rec3=3050;
    Rec4=240000;
    fa1=64./Re;
    fa2=19*Re.^-0.82;
    fb=(4.1*10^-16)*Re.^4;
    fc=0.351*Re.^-0.255;
    fd=0.118*Re.^-0.165;
    fatilde=fa1+(fa2-fa1)./(1+(Re/950).^-10).^0.5; 
    F1tilde=fatilde+(fb-fatilde)./(1+(Re./Rec).^-50).^0.5;
    F3tilde=F1tilde+(fc-F1tilde)./(1+(Re./Rec3).^-50).^0.5;
    lambda2=F3tilde+(fd-F3tilde)./(1+1./(Re/Rec4)).^0.5;
    ffr=lambda2;
end
ffrV(i)=ffr;
end

% figure(1)
% plot(ReV,ffrV,'b*')
% hold on
% grid on
% figure(2)
% semilogx(ReV,ffrV,'b*')
% grid on
% hold on