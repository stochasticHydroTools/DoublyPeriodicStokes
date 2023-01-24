%%

clear all
close all

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)
set(0,'defaultAxesLineWidth',3)

NBINS = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphalpha = 0.2;


EXP = load('./Exp_Data/Hist_Mat.mat');
b = EXP.b;
figure(1)
hplot(1) = plot(b,EXP.hm,'color',[255.0 180.0 40.0]./255.0);
hold all
pa = fill([b flip(b)], [EXP.hm+EXP.er fliplr(EXP.hm-EXP.er)],[255.0 180.0 40.0]./255.0);
set(pa,'facealpha',alphalpha,'edgecolor',[255.0 180.0 40.0]./255.0)
hold all

%%
phi = 0.4;
fps = 1; %1/10;
a = 1;
%close all
cols = parula(8);
% cols = [0.2161    0.7843    0.5923; 0.0689    0.6948    0.8394; 0.9970    0.7659    0.2199];
% cols = [[255.0 180.0 40.0]./255.0; [211.0 255.0 0.0]./255.0; [12.0 196.0 204.0]./255.0];
cols = [0.2422    0.1504    0.6603; 0.1786    0.5289    0.9682; 0.2161    0.7843    0.5923; 0.9778    0.0889    0.2667];
    
% Data = load('./MAT_FILES/LOWEST_FINAL_3_Torque_Lim.mat');

Data = load('./MAT_FILES/Spectral_torque_lim_Omega_FCM_pair_Case_Wall.mat');

k = 4;

figure(1)

hold all
hplot(3) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'color',cols(k,:));
hold all

Spect_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load('./MAT_FILES/Spectral_torque_lim_Omega_FCM_pair_Case_Wall_3xImages.mat');

k = 1;

figure(1)

hold all
hplot(4) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'--','color',cols(k,:));
hold all

Spect_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Data = load('./MAT_FILES/Spectral_torque_lim_Omega_FCM_pair_RPY_Wall.mat');
% 
% k = 2;
% 
% figure(1)
% 
% hold all
% hplot(3) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'color',cols(k,:));
% hold all
% 
% Spect_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data = load('./MAT_FILES/Spectral_torque_lim_Omega_COMPARE_OldC.mat');
% 
% k = 3;
% 
% figure(1)
% 
% hold all
% hplot(4) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'--','color',[1, 0.2, 0.55]);
% hold all
% 
% Spect_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load('./MAT_FILES/LOWEST_FINAL_3_Torque_Lim.mat');

k = 3;

figure(1)

hold all
hplot(2) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'-o','color',cols(k,:));
hold all

RPY_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);

xlim([0 80])
ylim([0 0.06])

legend(hplot,'Experiment','RPY','FCM, $$L=127 \, a$$',...
             'FCM, $$L=3 \times 127 \, a$$')
         
xlabel('$$U_x \ \mathrm{(\mu m / s)}$$')
ylabel('$$P(V)$$')
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

er_spec = sqrt(trapz(EXP.b,(Spect_interp-EXP.hm).^2))
er_rpy = sqrt(trapz(EXP.b,(RPY_interp-EXP.hm).^2))


%print('-dpng','-r300','/home/hat/GIT_Dirs/Papers/DPStokes/Figures/Rollers_Single_Wall.png')