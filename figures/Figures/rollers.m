clc;clear;close all
warning off
LOC = './data/rollers/';
format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.12;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5  ;
x_label_y = -0.13;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 4;
yl_in = 2.7;
dx_l_in = 0.57;dx_r_in = 0.1;
dy_b_in = 0.5 ;dy_t_in = 0.1;
xwidth = 1*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize   = 12 ;
fontsize       = 11 ;
small_fontsize = 10 ;
linewidth      = 1  ;
figbox_linewidth = 1;
ms = 4;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

% tikz parameters
image_ar = 9.17/3.71;
W_scaled = 0.7*xl;
W_in = W_scaled*xwidth;
H_in = W_in/image_ar;
x0 = (pos11(1)+pos11(3))*xwidth-W_in-0.795*figbox_linewidth*1/72;
y0 = (pos11(2)+pos11(4))*ywidth-1.1*H_in;
W = W_scaled*xwidth;

COLOR = {[255 180 40]/255,[0.2161 0.7843 0.5923],[0.9778 0.0889 0.2667],[0.2422 0.1504 0.6603]};

subplot('position',pos11)
hold on

alphalpha = 0.2;

EXP = load([LOC 'Hist_Mat.mat']);
b = EXP.b;
k = 1;
hand(k) = plot(b,EXP.hm,'color',COLOR{k},'linewidth',linewidth);
pa = fill([b flip(b)], [EXP.hm+EXP.er fliplr(EXP.hm-EXP.er)],COLOR{k});
set(pa,'facealpha',alphalpha,'edgecolor',COLOR{k})

Data = load([LOC 'LOWEST_FINAL_3_Torque_Lim.mat']);
RPY_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);
k = 2;
hand(k) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'-o','color',COLOR{k},'MarkerFaceColor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',ms,'linewidth',linewidth);

Data = load([LOC 'Spectral_torque_lim_Omega_FCM_pair_Case_Wall.mat']);
k = 3;
hand(k) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'color',COLOR{k},'linewidth',linewidth);

Data = load([LOC 'Spectral_torque_lim_Omega_FCM_pair_Case_Wall_3xImages.mat']);
Spect_interp = interp1(Data.SAVE_V_b,Data.SAVE_V_h,EXP.b);
k = 4;
hand(k) = plot(Data.SAVE_V_b,Data.SAVE_V_h,'-.','color',COLOR{k},'linewidth',linewidth);

% er_spec = sqrt(trapz(EXP.b,(Spect_interp-EXP.hm).^2));
% er_rpy = sqrt(trapz(EXP.b,(RPY_interp-EXP.hm).^2));

xlim([0,80])
ylim([0,0.057])
XLIM = xlim;
YLIM = ylim;
xlab = '$U_x\;(\mu\mathrm{m/s})$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\mathrm{P}\left(U_x\right)$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

x_dim = 20;
y_dim = 0.15;
clear text
text(x_dim,y_dim,'Monopole','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
x_dim = 20;
y_dim = 0.285;
clear text
text(x_dim,y_dim,'Dipole','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'Experiment','RPY','FCM, $L\!=\!127R_h$','FCM, $L\!=\!3\!\times\!127R_h$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.4*xl,dy_b+0.05*yl,0.05,0.05],'fontsize',small_fontsize);
leg.ItemTokenSize = [11,1];
leg.NumColumns = 2;
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[0 0.01 0.02 0.03 0.04 0.05 0.06],'xtick',[0 10 20 30 40 50 60 70 80]);
hold off

print('rollers','-dpdf','-painters')