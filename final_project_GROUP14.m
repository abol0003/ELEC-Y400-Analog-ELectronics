%% Analog Electronics final project
%%%%%
% Alexis BOLLENGIER 
% Sefa GÖNEN 
% GROUP 14
% 14/01/25
%%%%%
clc; 
close all; 
clear;

addpath(genpath('circuitDesign'));
addpath(genpath('models'));

load('models\UMC65_RVT.mat');

%% Initializations
designkitName		= 'umc65';
circuitTitle		= 'Analog Design - Project';
elementList.nmos	= {'Mn3','Mn4','Mn6'};
elementList.pmos	= {'Mp1','Mp2','Mp5','Mp7','Mp8'};

choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;

simulator			= 'spectre';
simulFile			= 0;
simulSkelFile		= 0;
spec				= [];
analog			= cirInit('analog',circuitTitle,'top',elementList,spec,choice,...
						designkitName,NRVT,PRVT,simulator,simulFile,simulSkelFile);
analog			= cirCheckInChoice(analog, choice);


%% Project: circuit
disp('                                                       ');
disp('  VDD          VDD                    VDD              ');
disp('   |             |                      |              ');
disp('  Mp8-+---------Mp7---------------------Mp5            ');
disp('   |--+          |                      |              ');
disp('   |          +--+--+         node 3->  +-----+---OUT  ');
disp('   |          |     |                   |     |        ');
disp('   |    IN1--Mp1   Mp2--IN2             |     |        ');
disp('   |          |     |                   |     |        ');
disp('   | node 1-> |     | <-node 2          |     Cl       ');
disp('   |          |--+  +------+-Cm---Rm----+     |        ');
disp(' Ibias        |  |  |      |    ↑       |     |        ');
disp('   |         Mn3-+-Mn4     |  node 4    |     |        ');
disp('   |          |     |      +-----------Mn6    |        ');
disp('   |          |     |                   |     |        ');
disp('  GND        GND   GND                 GND   GND       ');


%% AI: Implement your OpAmp according to your handcalculations
spec.VDD=1.1;
spec.FGBW=180e6;
spec.Cl = 10e-12;        
spec.Cm =spec.Cl/4 ;


%------------------Mn6----------------

Mn6.gm = 2*pi*spec.Cl*spec.FGBW;
Mn6.lg = 80e-9;
Mn6.vov = 0.0149;
Mn6.vds=spec.VDD/2;
Mn6.vsb=0; %triple-well tech assumptions
Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
Mn6.vgs=Mn6.vth+Mn6.vov;                   
Mn6.w = mosWidth('gm',Mn6.gm,Mn6); % width
Mn6 = mosNfingers(Mn6);
Mn6 = mosOpValues(Mn6);

%------------------Mp5 ------------------------------
Mp5.lg = 1e-6;
Mp5.ids = Mn6.ids; 
Mp5.vds=-spec.VDD/2;
Mp5.vov = -0.106;
Mp5.vsb=0; %triple-well technology assumptions
Mp5.vth= tableValueWref('vth', PRVT, Mp5.lg, 0, Mp5.vds, Mp5.vsb); 
Mp5.vgs =Mp5.vth+Mp5.vov;  %-0.325;
Mp5.w = mosWidth('ids', Mp5.ids, Mp5); 
Mp5   = mosNfingers(Mp5);
Mp5 = mosOpValues(Mp5);

%----------------Mp1, Mp2-------------------

Mp1.gm = 2*pi*spec.Cm*spec.FGBW;
Mp1.lg = 500e-9;
Mp1.vov = -0.120; 
Mp1.vds= -0.3695; %Impose value here first, Mp7 comes after
Mp1.vsb=0; %assumptions
Mp1.vth = tableValueWref('vth',PRVT,Mp1.lg,0,Mp1.vds,Mp1.vsb);
Mp1.vgs= Mp1.vth+Mp1.vov ; %Mp1.vgs=-0.3623; %-0.375 si pas d'interpolation
Mp1.w = mosWidth('gm',Mp1.gm,Mp1); % width
Mp1   = mosNfingers(Mp1);
Mp1 = mosOpValues(Mp1);


Mp2= cirElementCopy(Mp1, Mp2);


%------------------Mp7, Mp8 ------------------------------

Mp7.lg = 1e-6;
Mp7.ids = Mp1.ids*2; 
Mp7.vds=-0.3555;
Mp7.vsb = 0;      
Mp7.vgs=Mp5.vgs;
Mp7.w= mosWidth('ids', Mp7.ids, Mp7); 
Mp7= mosNfingers(Mp7); 
Mp7= mosOpValues(Mp7);

Mp8.lg = 1e-6;
Mp8.vgs = Mp5.vgs; 
Mp8.ids = Mp7.ids;
Mp8.vsb = 0;    %triple-well tech assumptions  
Mp8.vds =Mp8.vgs;
Mp8.w= mosWidth('ids', Mp8.ids, Mp8); 
Mp8= mosNfingers(Mp8); 
Mp8= mosOpValues(Mp8); 

%----------------Mn3, Mn4-----------
Mn3.lg = 0.5e-6;
Mn3.ids = Mp7.ids/2;
Mn3.vgs=Mn6.vgs;
Mn3.vds= Mn6.vgs;
Mn3.vsb=0; %triple-well tech assumptions
Mn3.w = mosWidth('ids',Mn3.ids,Mn3);
Mn3   = mosNfingers(Mn3);
Mn3 = mosOpValues(Mn3);

Mn4= cirElementCopy(Mn3, Mn4);


%% AI: Fill out the empty variables required to plot the transfer-function.
%  meaning of each variable see comment and
%  location of nodes see line 31 
spec.Rm = 5/Mn6.gm; 
%spec.Cm=0;

AvDC1 = Mp1.gm/(Mp1.gds+Mn3.gds) ;                           % DC gain 1st stage  19.2555;
AvDC2 = Mn6.gm/(Mn6.gds+Mp5.gds);                            % DC gain 2nd stage 16.4227;

%fig 8.3.3 of the chapter 8
Clprime=spec.Cl+Mn6.cdb+Mp5.cdb+Mp5.cgd;
Cn2=Mp2.cdb+Mn4.cdb+Mp2.cgd+Mn4.cgd+Mn6.cgs+Mn6.cgb;

C1    = Mp1.cgd +Mp1.cdb+ 2* Mn3.cgs+2*Mn3.cgb+Mn3.cgd;      % Capacitance on node 1
G1    = Mn3.gm+ Mn3.gds + Mp1.gds;                           % Admittance  on node 1

C2=Cn2+ spec.Cm * (1 + (Mn6.gm/(Mn6.gds + Mp5.gds)));       % Capacitance on node 2
G2    =  Mp2.gds + Mn4.gds;                                 % Admittance  on node 2

G3= Mn6.gm*spec.Cm;                                              % Admittance  on node 3
C3=spec.Cm*Clprime+spec.Cm*Cn2+Clprime*Cn2;                       % Capacitance on node 3

C4=Cn2;                                                           % Capacitance on node 4
G4=1/spec.Rm ;                                                    % Admittance  on node 4

% At very high frequencies,Cl and Cm behaves like a short circuit. 
% Explanation at end pg18 begin pg19 of chapter 8

%% AI: Set-up Rm, Cc and CL and calculate the zero required for the transfer-fct
spec.Rm =5/Mn6.gm;
spec.Gm=1/spec.Rm;
z2 = -(Mn6.gm/(4*spec.Cm)); % expression of zero due to the nulling resistor [rad/s]
fGBW = AvDC1*AvDC2*G2/(2*pi*C2);


%Optimal value
spec.Avtotal =  AvDC1*AvDC2; 
omega_pdom = spec.FGBW / spec.Avtotal; % Desired dominant pole
idealCm = (Mp2.gds + Mn4.gds)/(omega_pdom * (1+AvDC2));
idealRm = (1+ spec.Cl/idealCm)/Mn6.gm;
Avtotdb=20*log10(spec.Avtotal);

%% AI: Fill out the empty variables required for the performance summary
Vin_cm_min  =Mn3.vgs+abs(Mp1.vov)+Mp1.vgs ;   
Vin_cm_max  =spec.VDD+Mp7.vdsat+Mp1.vgs ; 
Vout_cm_min =Mn6.vdsat ;                         
Vout_cm_max =spec.VDD - abs(Mp5.vdsat) ;  
Pdiss       =spec.VDD*(Mp7.ids+Mn6.ids+ Mp8.ids) ;
Vcm=Vin_cm_min+((Vin_cm_max-Vin_cm_min)/2);

%% Sanity check (do not modify)

disp('======================================');
disp('=      Transistors in saturation     =');
disp('======================================');
if mosCheckSaturation(Mp1)
	fprintf('\nMp1:Success\n')
end
if mosCheckSaturation(Mp2)
	fprintf('Mp2:Success\n')
end
if mosCheckSaturation(Mn3)
	fprintf('Mn3:Success\n')
end
if mosCheckSaturation(Mn4)
	fprintf('Mn4:Success\n')
end
if mosCheckSaturation(Mp5)
	fprintf('Mp5:Success\n')
end
if mosCheckSaturation(Mn6)
	fprintf('Mn6:Success\n')
end
if mosCheckSaturation(Mp7)
	fprintf('Mp7:Success\n')
end
if mosCheckSaturation(Mp8)
	fprintf('Mp8:Success\n\n')
end


%% Summary of sizes and biasing points (do not modify)

disp('======================================');
disp('=    Sizes and operating points      =');
disp('======================================');
analog = cirElementsCheckOut(analog); % Update circuit file with 
% transistor sizes
mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% transistors in the circuit file
fprintf('IBIAS\t= %6.2fmA\nRm\t= %6.2f Ohm\nCm\t= %6.2fpF\n\n',Mp8.ids/1e-3,spec.Rm,spec.Cm/1e-12);

%% Performance summary (do not modify)

disp('======================================');
disp('=        Performance                 =');
disp('======================================');

fprintf('\nmetric        \t result\n');
fprintf('Vin,cm,min [mV] \t%.0f\n',Vin_cm_min/1e-3);
fprintf('Vin,cm,max [mV] \t%.0f\n',Vin_cm_max/1e-3);
fprintf('Vin,cm,    [mV] \t%.0f\n',Vcm/1e-3);
fprintf('Vout,cm,min [mV] \t%.0f\n',Vout_cm_min/1e-3);
fprintf('Vout,cm,max [mV] \t%.0f\n',Vout_cm_max/1e-3);
fprintf('Pdiss [mW]       \t%.1f\n\n',Pdiss/1e-3);

disp('======================================');
disp('=        Poles/Zeros                 =');
disp('======================================');

fprintf('\np1 = %.1f\tMHz\n', G1/(2*pi*C1*1e6));
fprintf('p2 = %.1f\tkHz\n', G2/(2*pi*C2*1e3));
fprintf('p3 = %.1f\tMHz\n', G3/(2*pi*C3*1e6));
fprintf('p4 = %.1f\tGHz\n', G4/(2*pi*C4*1e9));
fprintf('z1 = %.1f\tMHz\n', 2*G1/(2*pi*C1*1e6));
fprintf('z2 = %.1f\tMHz\n', z2/(-2*pi*1e6));

%% Ploting transfer function (do not modify)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
% transfer function
TF1   = AvDC1*AvDC2*((1+s*C1/(2*G1))*(1-s*(1/z2)))/ ...
                ((1+s*C1/G1)*(1+s*C2/G2)*(1+s*C3/G3)*(1+s*C4/G4));
[GM, PM, Wcg, Wcp] = margin(TF1);

freq = logspace(1,12,1e3);
figure(1)
bode(TF1,2*pi*freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
title(['Frequency response OTA with PM = ', num2str(PM), '°']);
hold all

spice_tbl = read_LTspice_data("final_project.txt");
[matlab_mag, matlab_phase] = bode(TF1, 2*pi*spice_tbl.freq);

figure
subplot(2, 1, 1);
semilogx(spice_tbl.freq*1e-6, 20*log10(squeeze(matlab_mag)));
hold on;
semilogx(spice_tbl.freq*1e-6, 20*log10(abs(spice_tbl.G1)));
grid on;
xlabel('Frequency (MHz)')
ylabel('Gain (dB)')
title('Full Bode Comparison')
subplot(2, 1, 2);
semilogx(spice_tbl.freq*1e-6, squeeze(matlab_phase));
hold on;
semilogx(spice_tbl.freq*1e-6, unwrap(angle(spice_tbl.G1))*180/pi);
grid on;
xlabel('Frequency (MHz)')
ylabel('Phase (°)')

fGBW = Mp2.gm / (2*pi*spec.Cm);
[~, minidx] = min(abs(2*fGBW - spice_tbl.freq));
roi = 1:(minidx+1);

figure
subplot(2, 1, 1);
semilogx(spice_tbl.freq(roi)*1e-6, 20*log10(squeeze(matlab_mag(:, :, roi))));
hold on;
semilogx(spice_tbl.freq(roi)*1e-6, 20*log10(abs(spice_tbl.G1(roi))));
grid on;
xlabel('Frequency (MHz)')
ylabel('Gain (dB)')
title('Zoomed-in Bode Comparison')
subplot(2, 1, 2);
semilogx(spice_tbl.freq(roi)*1e-6, squeeze(matlab_phase(:, :, roi)));
hold on;
semilogx(spice_tbl.freq(roi)*1e-6, unwrap(angle(spice_tbl.G1(roi)))*180/pi);
grid on;
xlabel('Frequency (MHz)')
ylabel('Phase (°)')

