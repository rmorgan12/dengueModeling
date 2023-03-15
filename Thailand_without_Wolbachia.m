% Using Ten Year Actual Temperature Data 
% Thailand

%
% Establishing Values
%

a=1;   % Start Day
b=4015; % End Day

t = a:b;   % Each Step

L=1.5;        % Ratio of Carrying Capacity over Human Pop.
gammaH=1/5;   % Progression from Exposed to Infected (Human)
sigma=.02;    % Recovery Rate
rhoN=1.25;    % Reproductive Rate (Vector)
rhoH=.00044;  % Reproductive Rate (Human)
eta=0.6228;   % Seasonal Amplitude
omega=20.61;  % Phase Shift
%phi = .001;  % Fatality Rate
muNO=1/14;    % Average Adult Death Rate Non-W
muH=.00001;   % Average Adult Death Rate Human Natural
reinfect = .08;
alpha=0.9;     % Maternal Transmission of Virus
h=0.4;  %Prey searching rate Assumed 0.1?0.65
alphaT=0.03; %Predator?s satiety rate Assumed 0.01?0.03
gammaT=0.15; %Predator?s conversation 
e=0.1; % Predator?s mortality rate 0.1
%
% Preallocation for Speed
%

tauN=zeros(1,length(t));muNA=zeros(1,length(t));
muN=zeros(1,length(t));bN=zeros(1,length(t));
c=zeros(1,length(t));Tmh=zeros(1,length(t));
Thm=zeros(1,length(t));

%
% Pulling Temperature Data
%

tempdata06 = rot90(xlsread('ThailandTemp.xlsx','A2:A366'));
tempdata07 = rot90(xlsread('ThailandTemp.xlsx','B2:B366'));
tempdata08 = rot90(xlsread('ThailandTemp.xlsx','C2:C366'));
tempdata09 = rot90(xlsread('ThailandTemp.xlsx','D2:D366'));
tempdata10 = rot90(xlsread('ThailandTemp.xlsx','E2:E366'));
tempdata11 = rot90(xlsread('ThailandTemp.xlsx','F2:F366'));
tempdata12 = rot90(xlsread('ThailandTemp.xlsx','G2:G366'));
tempdata13 = rot90(xlsread('ThailandTemp.xlsx','H2:H366'));
tempdata14 = rot90(xlsread('ThailandTemp.xlsx','I2:I366'));
tempdata15 = rot90(xlsread('ThailandTemp.xlsx','J2:J366'));
tempdata16 = rot90(xlsread('ThailandTemp.xlsx','K2:K366'));
tempAll = [tempdata06,tempdata07,tempdata08,tempdata09,tempdata10,tempdata11,tempdata12,tempdata13,tempdata14,tempdata15,tempdata16];
temp = tempAll(1:b);

%
% Temperature Dependant Factors
%

% Aquatic Vector Death Rate
for i=1:length(t)-1
    if (temp(i) >= 10.54) && (temp(i) <= 33.4)
        muNA(i) = .8692 - .159 * temp(i) + .01116 * temp(i).^2-.0003408 * temp(i).^3 + .000003809 * temp(i).^4;
    end
end

% Adult Vector Death Rate
for i=1:length(t)-1
    muN(i) = muNO*(1-eta*cos((2*pi*(t(i) + omega))/365)); 
end

% Biting Rate 
for i=1:length(t)-1
    if (temp(i)>=21) && (temp(i)<=32)
        bN(i) = .0943+.0043 * temp(i);
    else
        bN(i) = 0;
    end
end

% Transmission Chance Human to Vector
for i=1:length(t)
    if (temp(i) >= 12.4) && (temp(i) <= 26.1)
        Thm(i) = -.9037 + .0729 * temp(i);
    elseif (temp(i) < 12.4) || (temp(i)>32.5)
        Thm(i) = 0;
    elseif temp(i) > 26.1
        Thm(i)=1;
    end
end

 % Transmission Chance Vector to Human  %Is temp<=28 ? or temp<=32.5 ?
 for i=1:length(t)
     if (temp(i) >= 12.4) && (temp(i) <= 32.5)
         Tmh(i) = .001044 * temp(i) * (temp(i) - 12.286) * sqrt(32.461 - temp(i));
     else
         Tmh(i)=0;
     end
 end

% Vector Maturation Rate
for i=1:length(t)
    heavis = .57 - .43 * cos((2*pi*t(i))/365+(pi/4));
    if heavis<0
        heaviside=0;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    else
        heaviside=1;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    end
end

% Extrinsic Incubation Period
for i=1:length(t)
    c(i) = -.1393 + .008 * temp(i);
end

%
% Preallocation for Speed
%
FN=zeros(1,length(t));FH=zeros(1,length(t));SH=zeros(1,length(t));EH=zeros(1,length(t));
IH=zeros(1,length(t));RH=zeros(1,length(t));DH=zeros(1,length(t));AN=zeros(1,length(t));
SN=zeros(1,length(t));EN=zeros(1,length(t));IN=zeros(1,length(t));

%
% Initial Populations
%

% Human
sh=100000; % Suscetible
eh=0;      % Exposed
ih=1;   % Infected
rh=0;      % Recovered

% Non-W Vector
an=10000;  % Aquatic
sn=10000;  % Susceptible
en=0;      % Exposed
in=0;      % Infected

% Totals
th=sh+eh+ih+rh;
tn=sn+en+in;

% Normalize Initial Populations

SH(1) = sh/th;
EH(1) = eh/th;
IH(1) = ih/th;
RH(1) = rh/th;

AN(1) = an/tn;
SN(1) = sn/tn;
EN(1) = en/tn;
IN(1) = in/tn;

FN(1) = SN(1)+EN(1)+IN(1);
FH(1) = SH(1)+EH(1)+IH(1)+RH(1);

P(1) = 1;
% Differential Equations

for i = 1:length(t)-1
    FN(i) = SN(i) + EN(i) + IN(i);           % Non-W Pop. 
    FH(i) = SH(i) + EH(i) + IH(i) + RH(i);   % Human Pop.
    
    P(i+1) = P(i) + ((gammaT*h*(AN(i))/(1+alphaT*(AN(i)))))*P(i)-e*P(i);
    
    SH(i+1) = SH(i) + (- (bN(i) * Tmh(i) * L * IN(i) * SH(i)) + reinfect * RH(i)); 
    EH(i+1) = EH(i) + (bN(i) * Tmh(i) * L * IN(i) * SH(i) - gammaH * EH(i)); 
    IH(i+1) = IH(i) + (gammaH * EH(i) - sigma * IH(i));               
    RH(i+1) = RH(i) + (sigma * IH(i) - reinfect * RH(i));                                                
    
    AN(i+1) = AN(i) + (rhoN * (FN(i)/2) * (1 - AN(i))- (tauN(i) + muNA(i)) * AN(i))- P(i)*(h*AN(i))/(1+alphaT*AN(i)); 
    SN(i+1) = SN(i) + ((tauN(i) * AN(i)/2) - (bN(i) * Thm(i) * IH(i) + muN(i)) * SN(i)); 
    EN(i+1) = EN(i) + ((bN(i) * Thm(i) * IH(i)) * SN(i) - (c(i) + muN(i)) * EN(i)); 
    IN(i+1) = IN(i) + (c(i) * EN(i) - muN(i) * IN(i));                       
end
  
figure
subplot(4,2,[1,2])
plot(t,IH,'-r')
axis([1,b,0,.7])
title('Infected Human')

subplot(4,2,[3,4])
plot(t,FN,'-r')
axis([1,b,0,.7])
legend('total NW')
title('Mosquito Population')

%subplot(3,2,[3,4])
%plot(t,SH,'-b',t,EH,'-c',t,IH,'-r',t,RH,'-g')
%axis([1 b 0 1])
%legend('Susceptible', 'Exposed', 'Infected', 'Recovered')
%title('Human')

subplot(4,2,[7,8])
plot(t,temp)
axis([1 b 12 40])
title('Temperature')

subplot(4,2,[5,6])
plot(t,P,'-g')
axis([1 b 0 .7])
%legend('Toxochinites')
title('Predator')
