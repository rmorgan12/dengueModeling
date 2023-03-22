% Using Ten Year Actual Temperature Data
% Thailand
 
%
% Establishing Values
%
 
a=1;   % Start Day
b=4015; % End Day
 
t = a:b;   % Each Step
 
L=1.5;          % Ratio of Carrying Capacity over Human Pop.
gammaH=1/5;     % Progression from Exposed to Infected (Human)
gammaW = .1;    % Progression from Exposed to Infected (W-Vector)
sigma=.02;      % Recovery Rate
rhoN=1.25;      % Reproductive Rate (Vector)
rhoH=.00044;    % Reproductive Rate (Human)
eta=0.6228;     % Seasonal Amplitude
omega=20.61;    % Phase Shift
%phi = .001;    % Fatality Rate
muNO=1/14;      % Average Adult Death Rate Non-W
muH=.00001;     % Average Adult Death Rate Human Natural
muWA = 1/14;    % Wolbachia Aquatic Death Rate
reinfect=0.08; % Recovered to Susceptible
alpha=0.9;     % Maternal Transmission of Virus
h=0.65;  %Prey searching rate Assumed 0.1-0.65
alphaT=0.03; %Predator's satiety rate Assumed 0.01-0.03
gammaT=0.15; %Predator's conversation 
e=0.1; % Predator's mortality rate 0.1
 
%
% Preallocation for Speed
%
tauN=zeros(1,length(t));tauW=zeros(1,length(t));muNA=zeros(1,length(t));
muN=zeros(1,length(t));bN=zeros(1,length(t));
c=zeros(1,length(t));Tmh=zeros(1,length(t));
Thm=zeros(1,length(t));bW=zeros(1,length(t));
 
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
mean(temp)
std(temp)
%
% Temperature1	`	` Dependant Factors
%
 
% Aquatic Vector Death Rate
for i=1:length(t)-1
    if (temp(i) >= 10.54) && (temp(i) <= 33.4)
        muNA(i) = .8692 - .159 * temp(i) + .01116 * temp(i)^2 -.0003408 * temp(i)^3 + .000003809 * temp(i)^4;
    end
end
 
% Adult Vector Death Rate
for i=1:length(t)-1
    muN(i) = muNO * (1-eta*cos((2*pi*(t(i) + omega))/365)); 
end
 
% Biting Rate 
for i=1:length(t)-1
    if (temp(i)>=21) && (temp(i)<=32)
        bN(i) = .0943 +.0043 * temp(i);
        bW(i) = .95 * (.0943+.0043 * temp(i));
    else
        bN(i) = 0;
        bW(i) = 0;
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
        heaviside = 0;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    else
        heaviside = 1;        
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
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
SN=zeros(1,length(t));EN=zeros(1,length(t));IN=zeros(1,length(t));AW=zeros(1,length(t));
SW=zeros(1,length(t));EW=zeros(1,length(t));IW=zeros(1,length(t));FW=zeros(1,length(t));
 
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
 
% W-Vector
aw=10000;
sw=10000;
ew=0;
iw=0;
 
% Totals
th=sh+eh+ih+rh;
tn=sn+en+in;
tw=sw+ew+iw;
 
% Normalize Initial Populations
 
SH(1) = sh/th;
EH(1) = eh/th;
IH(1) = ih/th;
RH(1) = rh/th;
 
AN(1) = an/tn;
SN(1) = sn/tn;
EN(1) = en/tn;
IN(1) = in/tn;
 
AW(1) = aw/tw;
SW(1) = sw/tw;
EW(1) = ew/tw;
IW(1) = iw/tw;
 
FN(1) = SN(1)+EN(1)+IN(1);
FW(1) = SW(1)+EW(1)+IW(1);
FH(1) = SH(1)+EH(1)+IH(1)+RH(1);

P(1) = 1;

 
%
% Differential Equations
%
count = 0;
for i = 1:length(t)-1
    FN(i) = SN(i) + EN(i) + IN(i);          % Non-W Pop. 
    FW(i) = SW(i) + EW(i) + IW(i);          % W Pop. 
    FH(i) = SH(i) + EH(i) + IH(i) + RH(i);  % Human Pop.
    
    P(i+1) = P(i) + ((gammaT*h*(AN(i)+AW(i))/(1+alphaT*(AN(i)+AW(i)))))*P(i)-e*P(i);
     
    SH(i+1) = SH(i) + [-(bN(i) * Tmh(i) * L * IN(i) * SH(i)) - (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) + reinfect * RH(i)]; 
    EH(i+1) = EH(i) + [(bN(i) * Tmh(i) * L * IN(i) * SH(i)) + (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) - gammaH * EH(i)]; 
    IH(i+1) = IH(i) + [gammaH * EH(i) - sigma * IH(i)];               
    RH(i+1) = RH(i) + [sigma * IH(i) - reinfect * RH(i)];   
     
    AN(i+1) = AN(i) + [rhoN * (FN(i)^2)/(2*(FN(i) + FW(i))) * (1 - (AN(i) + AW(i)))- (tauN(i) + muNA(i)) * AN(i) - P(i)*(gammaT*h*AN(i))/(1+alphaT*(AN(i)+AW(i)))]; 
    SN(i+1) = SN(i) + [(tauN(i) * AN(i)/2) + (1 - alpha) * tauW(i) * AW(i)/2 - (bN(i) * FN(i) * IH(i) + muN(i)) * SN(i)]; 
    EN(i+1) = EN(i) + [(bN(i) * Thm(i) * IH(i) * SN(i)) - (c(i) + muN(i)) * EN(i)]; 
    IN(i+1) = IN(i) + [c(i) * EN(i) - muN(i) * IN(i)];    
     
    AW(i+1) = AW(i) + [.95 * rhoN * (FN(i)/2) * (1 - (AN(i) + AW(i)))- (tauW(i) + muWA) * AW(i) - (P(i)*(gammaT*h*AW(i))/(1+alphaT*(AN(i)+AW(i))))];
    SW(i+1) = SW(i) + [tauW(i) * alpha * (AW(i)/2) - (bW(i) * Thm(i) * IH(i) + 1.1 * muN(i)) * SW(i)];
    EW(i+1) = EW(i) + [(bW(i) * Thm(i) * IH(i) * SW(i)) - (gammaW + 1.1 * muN(i)) * EW(i)];
    IW(i+1) = IW(i) + [gammaW * EW(i) - 1.1 * muN(i) * IW(i)];

%   if count > 14
%       P(1+i) = P(i)+.25;
%       count = 0;
%   end
if i > 7
        TEMP = [temp(i-7),temp(i-6),temp(i-5),temp(i-4),temp(i-3),temp(i-2),temp(i-1)];
        if mean(TEMP) < 24 && count > 7
            P(i+1) = .25+P(i);
            count = 0;
         elseif P(i) < 0.1 && count > 7
             P(1+i) = .25 + P(i);
             count = 0;
        elseif count > 10
            P(i+1) = .25+P(i);
            count = 0;
        else
            count = count + 1;
        end
    end

 end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Using Ten Year Actual Temperature Data
% Thailand
 
%
% Establishing Values
%
 
a=1;   % Start Day
b=4015; % End Day
 
t = a:b;   % Each Step
 
L=1.5;          % Ratio of Carrying Capacity over Human Pop.
gammaH=1/5;     % Progression from Exposed to Infected (Human)
gammaW = .1;    % Progression from Exposed to Infected (W-Vector)
sigma=.02;      % Recovery Rate
rhoN=1.25;      % Reproductive Rate (Vector)
rhoH=.00044;    % Reproductive Rate (Human)
eta=0.6228;     % Seasonal Amplitude
omega=20.61;    % Phase Shift
%phi = .001;    % Fatality Rate
muNO=1/14;      % Average Adult Death Rate Non-W
muH=.00001;     % Average Adult Death Rate Human Natural
muWA = 1/14;    % Wolbachia Aquatic Death Rate
reinfect=0.08; % Recovered to Susceptible
alpha=0.9;     % Maternal Transmission of Virus
h=0.65;  %Prey searching rate Assumed 0.1-0.65
alphaT=0.03; %Predator's satiety rate Assumed 0.01-0.03
gammaT=0.15; %Predator's conversation 
e=0.1; % Predator's mortality rate 0.1
 
%
% Preallocation for Speed
%
tauN=zeros(1,length(t));tauW=zeros(1,length(t));muNA=zeros(1,length(t));
muN=zeros(1,length(t));bN=zeros(1,length(t));
c=zeros(1,length(t));Tmh=zeros(1,length(t));
Thm=zeros(1,length(t));bW=zeros(1,length(t));
 
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
% Temperature1	`	` Dependant Factors
%
 
% Aquatic Vector Death Rate
for i=1:length(t)-1
    if (temp(i) >= 10.54) && (temp(i) <= 33.4)
        muNA(i) = .8692 - .159 * temp(i) + .01116 * temp(i)^2 -.0003408 * temp(i)^3 + .000003809 * temp(i)^4;
    end
end
 
% Adult Vector Death Rate
for i=1:length(t)-1
    muN(i) = muNO * (1-eta*cos((2*pi*(t(i) + omega))/365)); 
end
 
% Biting Rate 
for i=1:length(t)-1
    if (temp(i)>=21) && (temp(i)<=32)
        bN(i) = .0943 +.0043 * temp(i);
        bW(i) = .95 * (.0943+.0043 * temp(i));
    else
        bN(i) = 0;
        bW(i) = 0;
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
        heaviside = 0;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    else
        heaviside = 1;        
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
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
IXH=zeros(1,length(t));RH=zeros(1,length(t));DH=zeros(1,length(t));AN=zeros(1,length(t));
SN=zeros(1,length(t));EN=zeros(1,length(t));IN=zeros(1,length(t));AW=zeros(1,length(t));
SW=zeros(1,length(t));EW=zeros(1,length(t));IW=zeros(1,length(t));FW=zeros(1,length(t));
 
%
% Initial Populations
%
 
% Human
sh=100000; % Suscetible
eh=0;      % Exposed
ixh=1;   % Infected
rh=0;      % Recovered
 
% Non-W Vector
an=10000;  % Aquatic
sn=10000;  % Susceptible
en=0;      % Exposed
in=0;      % Infected
 
% W-Vector
aw=10000;
sw=10000;
ew=0;
iw=0;
 
% Totals
th=sh+eh+ixh+rh;
tn=sn+en+in;
tw=sw+ew+iw;
 
% Normalize Initial Populations
 
SH(1) = sh/th;
EH(1) = eh/th;
IXH(1) = ixh/th;
RH(1) = rh/th;
 
AN(1) = an/tn;
SN(1) = sn/tn;
EN(1) = en/tn;
IN(1) = in/tn;
 
AW(1) = aw/tw;
SW(1) = sw/tw;
EW(1) = ew/tw;
IW(1) = iw/tw;
 
FN(1) = SN(1)+EN(1)+IN(1);
FW(1) = SW(1)+EW(1)+IW(1);
FH(1) = SH(1)+EH(1)+IXH(1)+RH(1);

X(1) = 1;

 
%
% Differential Equations
%
count = 0;
for i = 1:length(t)-1
    FN(i) = SN(i) + EN(i) + IN(i);          % Non-W Pop. 
    FW(i) = SW(i) + EW(i) + IW(i);          % W Pop. 
    FH(i) = SH(i) + EH(i) + IXH(i) + RH(i);  % Human Pop.
    
    X(i+1) = X(i) + ((gammaT*h*(AN(i)+AW(i))/(1+alphaT*(AN(i)+AW(i)))))*X(i)-e*X(i);
     
    SH(i+1) = SH(i) + [-(bN(i) * Tmh(i) * L * IN(i) * SH(i)) - (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) + reinfect * RH(i)]; 
    EH(i+1) = EH(i) + [(bN(i) * Tmh(i) * L * IN(i) * SH(i)) + (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) - gammaH * EH(i)]; 
    IXH(i+1) = IXH(i) + [gammaH * EH(i) - sigma * IXH(i)];               
    RH(i+1) = RH(i) + [sigma * IXH(i) - reinfect * RH(i)];   
     
    AN(i+1) = AN(i) + [rhoN * (FN(i)^2)/(2*(FN(i) + FW(i))) * (1 - (AN(i) + AW(i)))- (tauN(i) + muNA(i)) * AN(i) - X(i)*(gammaT*h*AN(i))/(1+alphaT*(AN(i)+AW(i)))]; 
    SN(i+1) = SN(i) + [(tauN(i) * AN(i)/2) + (1 - alpha) * tauW(i) * AW(i)/2 - (bN(i) * FN(i) * IXH(i) + muN(i)) * SN(i)]; 
    EN(i+1) = EN(i) + [(bN(i) * Thm(i) * IXH(i) * SN(i)) - (c(i) + muN(i)) * EN(i)]; 
    IN(i+1) = IN(i) + [c(i) * EN(i) - muN(i) * IN(i)];    
     
    AW(i+1) = AW(i) + [.95 * rhoN * (FN(i)/2) * (1 - (AN(i) + AW(i)))- (tauW(i) + muWA) * AW(i) - (X(i)*(gammaT*h*AW(i))/(1+alphaT*(AN(i)+AW(i))))];
    SW(i+1) = SW(i) + [tauW(i) * alpha * (AW(i)/2) - (bW(i) * Thm(i) * IXH(i) + 1.1 * muN(i)) * SW(i)];
    EW(i+1) = EW(i) + [(bW(i) * Thm(i) * IXH(i) * SW(i)) - (gammaW + 1.1 * muN(i)) * EW(i)];
    IW(i+1) = IW(i) + [gammaW * EW(i) - 1.1 * muN(i) * IW(i)];
    
%     if count > 14
%         X(i+1) = X(i)+.2;
%         count = 0;
%     end
    if i > 7
        TEMP = [temp(i-7),temp(i-6),temp(i-5),temp(i-4),temp(i-3),temp(i-2),temp(i-1)];
        if mean(TEMP) < 24 && count > 7
            X(i+1) = .2+X(i);
            count = 0;
         elseif X(i) < 0.1 && count > 7
             X(1+i) = .2 + X(i);
             count = 0;
        elseif count > 10
            X(i+1) = .2+X(i);
            count = 0;
        else
            count = count + 1;
        end
    end
end
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Ten Year Actual Temperature Data
% Thailand
 
%
% Establishing Values
%
 
a=1;   % Start Day
b=4015; % End Day
 
t = a:b;   % Each Step
 
L=1.5;          % Ratio of Carrying Capacity over Human Pop.
gammaH=1/5;     % Progression from Exposed to Infected (Human)
gammaW = .1;    % Progression from Exposed to Infected (W-Vector)
sigma=.02;      % Recovery Rate
rhoN=1.25;      % Reproductive Rate (Vector)
rhoH=.00044;    % Reproductive Rate (Human)
eta=0.6228;     % Seasonal Amplitude
omega=20.61;    % Phase Shift
%phi = .001;    % Fatality Rate
muNO=1/14;      % Average Adult Death Rate Non-W
muH=.00001;     % Average Adult Death Rate Human Natural
muWA = 1/14;    % Wolbachia Aquatic Death Rate
reinfect=0.08; % Recovered to Susceptible
alpha=0.9;     % Maternal Transmission of Virus
h=0.65;  %Prey searching rate Assumed 0.1-0.65
alphaT=0.03; %Predator's satiety rate Assumed 0.01-0.03
gammaT=0.15; %Predator's conversation 
e=0.1; % Predator's mortality rate 0.1
 
%
% Preallocation for Speed
%
tauN=zeros(1,length(t));tauW=zeros(1,length(t));muNA=zeros(1,length(t));
muN=zeros(1,length(t));bN=zeros(1,length(t));
c=zeros(1,length(t));Tmh=zeros(1,length(t));
Thm=zeros(1,length(t));bW=zeros(1,length(t));
 
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
% Temperature1	`	` Dependant Factors
%
 
% Aquatic Vector Death Rate
for i=1:length(t)-1
    if (temp(i) >= 10.54) && (temp(i) <= 33.4)
        muNA(i) = .8692 - .159 * temp(i) + .01116 * temp(i)^2 -.0003408 * temp(i)^3 + .000003809 * temp(i)^4;
    end
end
 
% Adult Vector Death Rate
for i=1:length(t)-1
    muN(i) = muNO * (1-eta*cos((2*pi*(t(i) + omega))/365)); 
end
 
% Biting Rate 
for i=1:length(t)-1
    if (temp(i)>=21) && (temp(i)<=32)
        bN(i) = .0943 +.0043 * temp(i);
        bW(i) = .95 * (.0943+.0043 * temp(i));
    else
        bN(i) = 0;
        bW(i) = 0;
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
        heaviside = 0;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    else
        heaviside = 1;        
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
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
IYH=zeros(1,length(t));RH=zeros(1,length(t));DH=zeros(1,length(t));AN=zeros(1,length(t));
SN=zeros(1,length(t));EN=zeros(1,length(t));IN=zeros(1,length(t));AW=zeros(1,length(t));
SW=zeros(1,length(t));EW=zeros(1,length(t));IW=zeros(1,length(t));FW=zeros(1,length(t));
 
%
% Initial Populations
%
 
% Human
sh=100000; % Suscetible
eh=0;      % Exposed
iyh=1;   % Infected
rh=0;      % Recovered
 
% Non-W Vector
an=10000;  % Aquatic
sn=10000;  % Susceptible
en=0;      % Exposed
in=0;      % Infected
 
% W-Vector
aw=10000;
sw=10000;
ew=0;
iw=0;
 
% Totals
th=sh+eh+iyh+rh;
tn=sn+en+in;
tw=sw+ew+iw;
 
% Normalize Initial Populations
 
SH(1) = sh/th;
EH(1) = eh/th;
IYH(1) = iyh/th;
RH(1) = rh/th;
 
AN(1) = an/tn;
SN(1) = sn/tn;
EN(1) = en/tn;
IN(1) = in/tn;
 
AW(1) = aw/tw;
SW(1) = sw/tw;
EW(1) = ew/tw;
IW(1) = iw/tw;
 
FN(1) = SN(1)+EN(1)+IN(1);
FW(1) = SW(1)+EW(1)+IW(1);
FH(1) = SH(1)+EH(1)+IYH(1)+RH(1);

Y(1) = 1;

 
%
% Differential Equations
%
count = 0;
for i = 1:length(t)-1
    FN(i) = SN(i) + EN(i) + IN(i);          % Non-W Pop. 
    FW(i) = SW(i) + EW(i) + IW(i);          % W Pop. 
    FH(i) = SH(i) + EH(i) + IYH(i) + RH(i);  % Human Pop.
    
    Y(i+1) = Y(i) + ((gammaT*h*(AN(i)+AW(i))/(1+alphaT*(AN(i)+AW(i)))))*Y(i)-e*Y(i);
     
    SH(i+1) = SH(i) + [-(bN(i) * Tmh(i) * L * IN(i) * SH(i)) - (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) + reinfect * RH(i)]; 
    EH(i+1) = EH(i) + [(bN(i) * Tmh(i) * L * IN(i) * SH(i)) + (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) - gammaH * EH(i)]; 
    IYH(i+1) = IYH(i) + [gammaH * EH(i) - sigma * IYH(i)];               
    RH(i+1) = RH(i) + [sigma * IYH(i) - reinfect * RH(i)];   
     
    AN(i+1) = AN(i) + [rhoN * (FN(i)^2)/(2*(FN(i) + FW(i))) * (1 - (AN(i) + AW(i)))- (tauN(i) + muNA(i)) * AN(i) - Y(i)*(gammaT*h*AN(i))/(1+alphaT*(AN(i)+AW(i)))]; 
    SN(i+1) = SN(i) + [(tauN(i) * AN(i)/2) + (1 - alpha) * tauW(i) * AW(i)/2 - (bN(i) * FN(i) * IYH(i) + muN(i)) * SN(i)]; 
    EN(i+1) = EN(i) + [(bN(i) * Thm(i) * IYH(i) * SN(i)) - (c(i) + muN(i)) * EN(i)]; 
    IN(i+1) = IN(i) + [c(i) * EN(i) - muN(i) * IN(i)];    
     
    AW(i+1) = AW(i) + [.95 * rhoN * (FN(i)/2) * (1 - (AN(i) + AW(i)))- (tauW(i) + muWA) * AW(i) - (Y(i)*(gammaT*h*AW(i))/(1+alphaT*(AN(i)+AW(i))))];
    SW(i+1) = SW(i) + [tauW(i) * alpha * (AW(i)/2) - (bW(i) * Thm(i) * IYH(i) + 1.1 * muN(i)) * SW(i)];
    EW(i+1) = EW(i) + [(bW(i) * Thm(i) * IYH(i) * SW(i)) - (gammaW + 1.1 * muN(i)) * EW(i)];
    IW(i+1) = IW(i) + [gammaW * EW(i) - 1.1 * muN(i) * IW(i)];
    
%     if count > 14
%         Y(1+i) = Y(i) +.15;
%         count = 0;
%     end
    if i > 7
        TEMP = [temp(i-7),temp(i-6),temp(i-5),temp(i-4),temp(i-3),temp(i-2),temp(i-1)];
        if mean(TEMP) < 24 && count > 7
            Y(i+1) = .15+Y(i);
            count = 0;
         elseif Y(i) < 0.1 && count > 7
             Y(1+i) = .15 + Y(i);
             count = 0;
        elseif count > 10
            Y(i+1) = .15+Y(i);
            count = 0;
        else
            count = count + 1;
        end
    end
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using Ten Year Actual Temperature Data
% Thailand
 
%
% Establishing Values
%
 
a=1;   % Start Day
b=4015; % End Day
 
t = a:b;   % Each Step
 
L=1.5;          % Ratio of Carrying Capacity over Human Pop.
gammaH=1/5;     % Progression from Exposed to Infected (Human)
gammaW = .1;    % Progression from Exposed to Infected (W-Vector)
sigma=.02;      % Recovery Rate
rhoN=1.25;      % Reproductive Rate (Vector)
rhoH=.00044;    % Reproductive Rate (Human)
eta=0.6228;     % Seasonal Amplitude
omega=20.61;    % Phase Shift
%phi = .001;    % Fatality Rate
muNO=1/14;      % Average Adult Death Rate Non-W
muH=.00001;     % Average Adult Death Rate Human Natural
muWA = 1/14;    % Wolbachia Aquatic Death Rate
reinfect=0.08; % Recovered to Susceptible
alpha=0.9;     % Maternal Transmission of Virus
h=0.65;  %Prey searching rate Assumed 0.1-0.65
alphaT=0.03; %Predator's satiety rate Assumed 0.01-0.03
gammaT=0.15; %Predator's conversation 
e=0.1; % Predator's mortality rate 0.1
 
%
% Preallocation for Speed
%
tauN=zeros(1,length(t));tauW=zeros(1,length(t));muNA=zeros(1,length(t));
muN=zeros(1,length(t));bN=zeros(1,length(t));
c=zeros(1,length(t));Tmh=zeros(1,length(t));
Thm=zeros(1,length(t));bW=zeros(1,length(t));
 
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
% Temperature1	`	` Dependant Factors
%
 
% Aquatic Vector Death Rate
for i=1:length(t)-1
    if (temp(i) >= 10.54) && (temp(i) <= 33.4)
        muNA(i) = .8692 - .159 * temp(i) + .01116 * temp(i)^2 -.0003408 * temp(i)^3 + .000003809 * temp(i)^4;
    end
end
 
% Adult Vector Death Rate
for i=1:length(t)-1
    muN(i) = muNO * (1-eta*cos((2*pi*(t(i) + omega))/365)); 
end
 
% Biting Rate 
for i=1:length(t)-1
    if (temp(i)>=21) && (temp(i)<=32)
        bN(i) = .0943 +.0043 * temp(i);
        bW(i) = .95 * (.0943+.0043 * temp(i));
    else
        bN(i) = 0;
        bW(i) = 0;
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
        heaviside = 0;
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
    else
        heaviside = 1;        
        tauN(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
        tauW(i) = (.00483 * temp(i) -.00796)*(.57-.43*cos((2*pi*t(i))/365+(pi/4)))*heaviside;
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
IZH=zeros(1,length(t));RH=zeros(1,length(t));DH=zeros(1,length(t));AN=zeros(1,length(t));
SN=zeros(1,length(t));EN=zeros(1,length(t));IN=zeros(1,length(t));AW=zeros(1,length(t));
SW=zeros(1,length(t));EW=zeros(1,length(t));IW=zeros(1,length(t));FW=zeros(1,length(t));
 
%
% Initial Populations
%
 
% Human
sh=100000; % Suscetible
eh=0;      % Exposed
izh=1;   % Infected
rh=0;      % Recovered
 
% Non-W Vector
an=10000;  % Aquatic
sn=10000;  % Susceptible
en=0;      % Exposed
in=0;      % Infected
 
% W-Vector
aw=10000;
sw=10000;
ew=0;
iw=0;
 
% Totals
th=sh+eh+izh+rh;
tn=sn+en+in;
tw=sw+ew+iw;
 
% Normalize Initial Populations
 
SH(1) = sh/th;
EH(1) = eh/th;
IZH(1) = izh/th;
RH(1) = rh/th;
 
AN(1) = an/tn;
SN(1) = sn/tn;
EN(1) = en/tn;
IN(1) = in/tn;
 
AW(1) = aw/tw;
SW(1) = sw/tw;
EW(1) = ew/tw;
IW(1) = iw/tw;
 
FN(1) = SN(1)+EN(1)+IN(1);
FW(1) = SW(1)+EW(1)+IW(1);
FH(1) = SH(1)+EH(1)+IZH(1)+RH(1);

Z(1) = 1;

 
%
% Differential Equations
%
count = 0;
for i = 1:length(t)-1
    FN(i) = SN(i) + EN(i) + IN(i);          % Non-W Pop. 
    FW(i) = SW(i) + EW(i) + IW(i);          % W Pop. 
    FH(i) = SH(i) + EH(i) + IZH(i) + RH(i);  % Human Pop.
    
    Z(i+1) = Z(i) + ((gammaT*h*(AN(i)+AW(i))/(1+alphaT*(AN(i)+AW(i)))))*Z(i)-e*Z(i);
     
    SH(i+1) = SH(i) + [-(bN(i) * Tmh(i) * L * IN(i) * SH(i)) - (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) + reinfect * RH(i)]; 
    EH(i+1) = EH(i) + [(bN(i) * Tmh(i) * L * IN(i) * SH(i)) + (bW(i) * .5 * Tmh(i) * L * IW(i) * SH(i)) - gammaH * EH(i)]; 
    IZH(i+1) = IZH(i) + [gammaH * EH(i) - sigma * IZH(i)];               
    RH(i+1) = RH(i) + [sigma * IZH(i) - reinfect * RH(i)];   
     
    AN(i+1) = AN(i) + [rhoN * (FN(i)^2)/(2*(FN(i) + FW(i))) * (1 - (AN(i) + AW(i)))- (tauN(i) + muNA(i)) * AN(i) - Z(i)*(gammaT*h*AN(i))/(1+alphaT*(AN(i)+AW(i)))]; 
    SN(i+1) = SN(i) + [(tauN(i) * AN(i)/2) + (1 - alpha) * tauW(i) * AW(i)/2 - (bN(i) * FN(i) * IZH(i) + muN(i)) * SN(i)]; 
    EN(i+1) = EN(i) + [(bN(i) * Thm(i) * IZH(i) * SN(i)) - (c(i) + muN(i)) * EN(i)]; 
    IN(i+1) = IN(i) + [c(i) * EN(i) - muN(i) * IN(i)];    
     
    AW(i+1) = AW(i) + [.95 * rhoN * (FN(i)/2) * (1 - (AN(i) + AW(i)))- (tauW(i) + muWA) * AW(i) - (Z(i)*(gammaT*h*AW(i))/(1+alphaT*(AN(i)+AW(i))))];
    SW(i+1) = SW(i) + [tauW(i) * alpha * (AW(i)/2) - (bW(i) * Thm(i) * IZH(i) + 1.1 * muN(i)) * SW(i)];
    EW(i+1) = EW(i) + [(bW(i) * Thm(i) * IZH(i) * SW(i)) - (gammaW + 1.1 * muN(i)) * EW(i)];
    IW(i+1) = IW(i) + [gammaW * EW(i) - 1.1 * muN(i) * IW(i)];
    
%     if count > 14
%         Z(i+1) = Z(i) + .1;
%         count = 0;
%     end
    if i > 7
        TEMP = [temp(i-7),temp(i-6),temp(i-5),temp(i-4),temp(i-3),temp(i-2),temp(i-1)];
        if mean(TEMP) < 24 && count > 7
            Z(i+1) = .1+Z(i);
            count = 0;
         elseif Z(i) < 0.1 && count > 7
             Z(1+i) = .1 + Z(i);
             count = 0;
        elseif count > 10
            Z(i+1) = .1+Z(i);
            count = 0;
        else
            count = count + 1;
        end
    end
end
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
subplot(2,2,[1,2])
plot(t,IH,'-b',t,IXH,'-r',t,IYH,'-g',t,IZH,'-m')
axis([1 b 0 0.7])
xlabel('Days')
legend({'25% Predator','20% Predator','15% Predator','10% Predator'},'FontSize',14)
title('Infected Human', 'fontsize', 18)

subplot(2,2,[3,4])
plot(t,AN,'-r',t,AW,'-b', t, P, '-g')
axis([1 b 0 1])
legend({'total AN','total AM', 'total P'},'FontSize',14)
xlabel('Days')
title('Aquatic Mosquito Population', 'fontsize', 18)
 
% subplot(3,2,[5,6])
% plot(t,temp,'-b')
% axis([1 b 12 40])
% legend({'Temperature'},'FontSize',14)
% title('Temp')
