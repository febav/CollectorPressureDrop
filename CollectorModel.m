%%
%     GNU General Public License
%   
%     CollectorPressureDrop is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     CollectorPressureDrop is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%     You should have received a copy of the GNU General Public License
%     along with CollectorPressureDrop. If not, see <http://www.gnu.org/licenses/>.
% 
%     Author:
%     Federico Bava
%     PhD student at the Technical University of Denmark (DTU)
%     febav@byg.dtu.dk
%  
%     Copyright (c) 2015 Federico Bava 
%     Released under GNU GPL v3

%%   CollectorModel.m  solve for the Flow Rates in two Parallel Pipes 
%
%   This file computes the flow distribution and pressure drop in a simple
%   parallel piping system, iteratively solving a set of nonlinear
%   equations, where the coefficient matrix is function of the solution vector.
%   Given an initial guess of the flow distribution Qold, the coefficient
%   matrix AA can be determined and the solution of a linear system,
%   AA(x)*x = BB, can be obtained. This process is continued until the user
%   defined tolerance is reached or the number of iterations has exceeded.
%   The structure of the matrix problem is  like:
%
%   (M1new)   (   1         1         1         1         1      | Mtot )
%   (M2new)   ( Y1*M1old -Y2*M2old    0         0         0      |  0   )
%   (M3new) = (   0       Y2*M2old -Y3*M3old    0         0      |  0   )
%   (M4new)   (   0         0       Y3*M3old -Y4*M4old    0      |  0   )
%   (M5new)   (   0         0         0       Y4*M4old -Y5*M5old |  0   )
%
%   where   M is the mass flow rate
%           tot denotes the total flow rate across the collector
%           Y is pressure drop resistance coefficient
%           1,2,...5 denote the hydraulic path
%           new denotes the flow rate computed at the present iteration
%           old denotes the flow rate computed at the previous iteration

%% code
clear all
close all
clc

% set parallel parameters
Np=18;                    % number of hydraulic paths (horizontal pipes)

% set pipe parameters 
DCi=0.0329*ones(Np,1);  % [m] diameter of inlet manifold
DCo=0.0329*ones(Np,1);  % [m] diameter of outlet manifold
DR=0.0073*ones(Np,1);   % [m] diameter of fictitious row pipe
ACi=DCi.^2*pi/4;        % [m^2] flow area of main pipes supplying each row
ACo=DCo.^2*pi/4;        % [m^2] flow area of main pipes returning from each row
AR=DR.^2*pi/4;          % [m^2] flow area of horizontal pipe 
LC=[0.165;0.1215*ones(Np-1,1)]; % [m] distance between horizontal pipes
LR=5.8*ones(Np,1);      % [m] length of horizontal pipes
ep=0;                   % [m] surface roughness
frfac_correlation=1;    % 1=Blasius, 2=Colebrook, 3=Haaland. 4=Joseph&Jang

if Np~=length(DCo) && Np~=length(DCi) && Np~=length(DR) && Np~=length(LR)
    error('Np is different from specified number of rows!')
end

% initialization
fid = 1;            % input for fprintf function (1 = print on screen)
it=1;
itmax=40;           % max number of iterations
itend=itmax;
tolFlow=1e-5;       % tolerance for flow convergence
tolDp=1e-4;         % tolerance for pressure convergence
ReR=zeros(Np,1);    % [-] Reynolds number, R postscript refers to absorber pipe
fR=zeros(Np,1);     % [-] Darcy friction factor, // 
YR=zeros(Np,1);     % [1/kg.m] resist. coeff.,   //
ReCi=zeros(Np,1);   % C postscript refers to manifold segments between pipes
ReCo=zeros(Np,1);   % i/o refers to inlet/outlet manifold
fCi=zeros(Np,1);    % //
fCo=zeros(Np,1);    % //
YTinSi=zeros(Np,1); % Inlet=diverter press.drop coeff. (side), Dp=Y_CSi*M^2, M is the flow above the diverter
YTinSt=zeros(Np,1); % Inlet=diverter  //        // (straight),  //
ZTinSi=zeros(Np,1); % Inlet=diverter  //        // (side), Idelchick definition
ZTinSt=zeros(Np,1); % Inlet=diverter  //        // (straight),  //
YToutSi=zeros(Np,1); % Out=converter press.drop coeff. (side), Dp=Y_CSi*M^2, M is the flow downst. the converter
YToutSt=zeros(Np,1); % Out=converter  //        // (straight),  //
ZToutSi=zeros(Np,1); % Out=converter  //        // (side), Idelchick definition
ZToutSiL=zeros(Np,1);% Out=converter  //        // (side+laminar), Idelchick definition
ZToutSt=zeros(Np,1); % Out=converter  //        // (straight),  //
YCi=zeros(Np,1);    % YCi/o(1) refers to the segment of manifold supplying pipe1. DpC=YCi*MC^2;
YCo=zeros(Np,1);    % YCi/o(2) refers to the segment between absorber pipe 1 and 2, and so on
YTot=zeros(Np,1);   % T postscript refers to "total" ( R + sum(C_i) )

% fluid properties
x=40;             % [%] glycol content (=0 for water)
Tin=20;     % [degC] fluid (inlet) temperature 
rhoIn=densityGlyMixAndWat(x,Tin);  % [kg/m3] fluid density at inlet
nuIn=viscosityGlyMixAndWat(x,Tin); % [m^2/s] fluid kin. viscosity at inlet 
Vtot_m3h=2.4;         % [m3/h] total volume flow rate in input
Vtot=Vtot_m3h/3600;     % [m3/s] total volume rate in input
Mtot=Vtot*rhoIn;        % [kg/s] total mass rate in input

% set reasonable initial guess for flow rates
MoldR=Mtot/Np*ones(Np,1); % [kg/s] mass flow rate per row (uniform distribution assumed)
MoldC=Mtot*ones(Np,1);    % MoldC(1)=Mtot (other cells will be overwritten)
for jj=2:Np 
  MoldC(jj)=MoldC(jj-1)-MoldR(jj-1); % mass flow in the connecting pipes 
end
MnewC=MoldC; % creation of variable MnewC, so that the first cell is Mtot

% set right hand side (RHS) column vector
BB=[Mtot zeros(1,Np-1)]';    % RHS vector

while it<=itmax
    figure(1)
    set(gcf,'DefaultAxesColorOrder',jet(20));
    figure(1)
    plot(1:Np,MoldR/Mtot*100)
    xlabel('Row #')
    ylabel('% of total flow')
    grid on
    hold all
    
    VoldR=MoldR/rhoIn;  % [m3/s] volume flow rate in absorber pipes
    VoldC=MoldC/rhoIn;  % [m3/s] volume flow rate in manifold segments
  for jj=1:Np   
    ReR(jj)=(DR(jj)./(nuIn*AR(jj))).*abs(MoldR(jj)/rhoIn);  % Re in abs. pipes [-]
    fR(jj)=FrictionFactorFunc(ReR(jj),DR(jj),ep,frfac_correlation); % fr.f. in abs. pipes [-]
    YR(jj)=fR(jj).*LR(jj)./DR(jj)*0.5/rhoIn./(AR(jj).^2); % resistance coefficients [1/kg.m]
    ReCi(jj)=(DCi(jj)./(nuIn*ACi(jj))).*abs(MoldC(jj)/rhoIn);   % Re in manifold in  [-]
    ReCo(jj)=(DCo(jj)./(nuIn*ACo(jj))).*abs(MoldC(jj)/rhoIn);   % Re in manifold out [-]
    fCi(jj)=FrictionFactorFunc(ReCi(jj),DCi(jj),ep,frfac_correlation); % fr.f. in manifold in [-]
    fCo(jj)=FrictionFactorFunc(ReCo(jj),DCo(jj),ep,frfac_correlation); % fr.f. in manifold out [-]
    [YTinSi(jj),ZTinSi(jj)]=TeeDivSide(VoldR(jj),VoldC(jj),AR(jj),ACi(jj),ACi(jj),rhoIn,ReCi(jj));  % diverter=inlet (side)
    [YTinSt(jj),ZTinSt(jj)]=TeeDivSt(VoldR(jj),VoldC(jj),AR(jj),ACi(jj),ACi(jj),rhoIn,ReCi(jj));  % diverter=inlet (straight)
    [YToutSi(jj),ZToutSi(jj),ZToutSiL(jj)]=TeeConvSide(VoldR(jj),VoldC(jj),AR(jj),ACi(jj),ACi(jj),rhoIn,ReCo(jj)); % converter=out (side)
    [YToutSt(jj),ZToutSt(jj)]=TeeConvSt(VoldR(jj),VoldC(jj),AR(jj),ACi(jj),ACi(jj),rhoIn,ReCo(jj),ZToutSiL(jj));   % converter=out (straight)
    YCi(jj)=fCi(jj).*LC(jj)./DCi(jj)*0.5/rhoIn./(ACi(jj).^2); % res.coef. in [1/kg.m]
    YCo(jj)=fCo(jj).*LC(jj)./DCo(jj)*0.5/rhoIn./(ACo(jj).^2); % res.coef. out [1/kg.m]   
    YTot(jj)=YR(jj);    % YT inizialized as YR, as it includes at least the row Dp
    for ii=1:jj
      YTot(jj)=YTot(jj)+...                             % row
          YCi(ii)*MoldC(ii)^2/MoldR(jj)^2+... % supply pipe
          YCo(ii)*MoldC(ii)^2/MoldR(jj)^2+... % return pipe
          (ii<jj)*YTinSt(ii)*MoldC(ii)^2/MoldR(jj)^2+...  % T-in (div) straight, normalized for MC
          (ii<jj)*YToutSt(ii)*MoldC(ii)^2/MoldR(jj)^2;  % T-out(con) straight, normalized for MC
    end
    YTot(jj)=YTot(jj)+YTinSi(jj)*MoldC(jj)^2/MoldR(jj)^2+...
                     +YToutSi(jj)*MoldC(jj)^2/MoldR(jj)^2; % add side passage of the inlet-tee
  end
  
  AA=diag(-YTot.*MoldR,0)+diag(YTot(1:end-1).*MoldR(1:end-1),-1);
  AA(1,:)=1;              % coefficent matrix
  MnewR = AA\BB;          % find solution vector
  for jj=2:Np
    MnewC(jj)=MnewC(jj-1)-MnewR(jj-1);
  end

  % check convergence
  Dp=[MnewR(1:end-1).*diag(AA,-1);MnewR(end)*(-AA(end,end))];
  emaxR=max(abs(MnewR-MoldR)./MnewR);
  emaxC=max(abs(MnewC-MoldC)./MnewC);
  emaxFlow=max(emaxC,emaxR);
  emaxDp=(max(Dp)-min(Dp))/mean(Dp);
  
  fprintf(1,'\nIteration = %3d      max error = %8.2e \n',it,emaxFlow);
  fprintf(1,'MnewR [kg/s]  MoldR [kg/s]  MnewC [kg/s]  MoldC [kg/s]\n');
  for j = 1:length(MoldR)
    fprintf(1,'%8.2g    %11.2g  %11.2g  %11.2g\n',MnewR(j),MoldR(j),MnewC(j),MoldC(j));
  end
 
  if emaxFlow<tolFlow && emaxDp<tolDp
    itend=it;
    it=itmax+1; 
  else 
    it=it+1;
    MoldR=(MnewR+MoldR)/2;
    for jj=2:Np
      MoldC(jj)=MoldC(jj-1)-MoldR(jj-1);
    end
  end
end

%   print max relative error and iteration count
fprintf(1,'\n Number of iterations to convergence = %3d\n',itend);
fprintf(1,' Max relative error at convergence =  %8.2e\n',emaxFlow);
if itend==itmax 
  warning(' WARNING  --  Hit max number of iterations!!! \n');
end

VR=MnewR/rhoIn;        % [m3/s] volume flow rate in different rows
vR=VR./AR;          % [m/s] average flow velocity in rows
VC=MnewC/rhoIn;          % [m3/s] volume flow rate in connecting pipes
vCi=VC./ACi;        % [m/s] average flow velocity conn. pipes (inlet)
vCo=VC./ACi;        % [m/s] average flow velocity conn. pipes (outlet)
DpR=(MnewR.^2).*YR;    % [Pa] pressure drop in row only
DpCi=(MnewC.^2).*YCi;  % [Pa] Pressure drop in conn. pipe segment
DpCo=(MnewC.^2).*YCo;  % [Pa] Pressure drop in conn. pipe segment
 
% print summary results
fprintf(fid,'\n ParallelFlow.m:  Summary Results \n\n');
fprintf(fid,'   Calculated Parameters for each Row: \n');
fprintf(fid,'Pipe Mass Rate Flow Rate Velocity  Reynolds Frict.Loss\n');
fprintf(fid,'      [kg/s]    [m^3/h]   [m/s]     [-]       [Pa]  \n');
for n=1:length(MnewR)
fprintf(fid,'%2i  %7.2f   %7.3f  %7.2f  %8.0f %9.0f \n', ...
            n,MnewR(n),VR(n)*3600,vR(n),ReR(n),DpR(n));
end
fprintf(fid,' \n');
fprintf(fid,'Fluid Properties: \n');
fprintf(fid,' density [kg/m^3]:          %10.3e \n',rhoIn);
fprintf(fid,' kinem. viscosity [m^2/s]:  %10.3e \n',nuIn);
fprintf(fid,' inlet temperature [degC]:  %10.1f \n',Tin);
fprintf(fid,' glycol content [%%]:        %10.1f \n',x);
fprintf(fid,' total flow rate [kg/s]:    %10.4f \n',Mtot);
fprintf(fid,' total flow rate [m3/h]:    %10.4f \n',Vtot*3600);
fprintf(fid,' Total Calculated Flow Rate [kg/s]:               %7.4f \n',sum(MnewR));
fprintf(fid,' Calculated DeltaP across Branched Section [kPa]: %7.2f \n\n',Dp(1)/1000);
