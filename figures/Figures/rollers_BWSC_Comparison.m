clc;clear;close all
warning off

format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.12;
y_label_y = 0.5 ;
y_label_x_inset = -0.17;
y_label_y_inset = 0.5 ;

x_label_x = 0.5  ;
x_label_y = -0.13;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.2;

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
linewidth      = 1.5;
figbox_linewidth = 1;
ms = 5;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;
fac = 2;
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

%%
Data = load('data/rollers/MAT_FILES/Spectral_torque_lim_Omega_FCM_pair_Case_Wall_3xImages_hist.mat');

Wall_b = Data.SAVE_V_b;
Wall_v = Data.SAVE_V_h;
Wall_Z = Data.SAVE_Z_b;
Wall_hZ = Data.SAVE_Z_h;
Wall_b_high = Data.SAVE_V_b_high;
Wall_V_high = Data.SAVE_V_h_high;
Wall_b_low = Data.SAVE_V_b_low;
Wall_V_low = Data.SAVE_V_h_low;

ref = 1;
fac = 4;
if ref
    [x_ref,y_ref] = SPLINE(Wall_b,Wall_v,fac);
    hand(1) = plot(x_ref,y_ref,'-','color',COLOR{1},'linewidth',linewidth);
    [x_ref,y_ref] = SPLINE(Wall_b_low,Wall_V_low,fac);
    hand(2) = plot(x_ref,y_ref,'--','color',COLOR{1},'linewidth',linewidth);
    [x_ref,y_ref] = SPLINE(Wall_b_high,Wall_V_high,fac);
    hand(3) = plot(x_ref,y_ref,':','color',COLOR{1},'linewidth',linewidth);
else
    hand(1) = plot(Wall_b,Wall_v,'-','color',COLOR{1},'linewidth',linewidth);
    hand(2) = plot(Wall_b_low,Wall_V_low,'--','color',COLOR{1},'linewidth',linewidth);
    hand(3) = plot(Wall_b_high,Wall_V_high,':','color',COLOR{1},'linewidth',linewidth);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load('data/rollers/MAT_FILES/Spectral_Channel_H6a.mat');

Chan6a_b = Data.SAVE_V_b;
Chan6a_v = Data.SAVE_V_h;
Chan6a_Z = Data.SAVE_Z_b;
Chan6a_hZ = Data.SAVE_Z_h;
Chan6a_b_high = Data.SAVE_V_b_high;
Chan6a_V_high = Data.SAVE_V_h_high;
Chan6a_b_low = Data.SAVE_V_b_low;
Chan6a_V_low = Data.SAVE_V_h_low;

if ref
    [x_ref,y_ref] = SPLINE(Chan6a_b,Chan6a_v,fac);
    hand(4) = plot(x_ref,y_ref,'-','color',COLOR{2},'linewidth',linewidth);
    [x_ref,y_ref] = SPLINE(Chan6a_b_low,Chan6a_V_low,fac);
    hand(5) = plot(x_ref,y_ref,'--','color',COLOR{2},'linewidth',linewidth);
    [x_ref,y_ref] = SPLINE(Chan6a_b_high,Chan6a_V_high,fac);
    hand(6) = plot(x_ref,y_ref,':','color',COLOR{2},'linewidth',linewidth);
else
    hand(4) = plot(Chan6a_b,Chan6a_v,'-','color',COLOR{2},'linewidth',linewidth);
    hand(5) = plot(Chan6a_b_low,Chan6a_V_low,'--','color',COLOR{2},'linewidth',linewidth);
    hand(6) = plot(Chan6a_b_high,Chan6a_V_high,':','color',COLOR{2},'linewidth',linewidth);
end

xlim([-5 60])
ylim([0,0.06])

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

LEGEND = {'BW','$h<2R_h$','$h>2R_h$','SC','$h<2R_h$','$h>2R_h$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.75*xl,dy_b+0.09*yl,0.05,0.05],'fontsize',small_fontsize);
leg.ItemTokenSize = [14,1];
leg.NumColumns = 2;
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[0 0.01 0.02 0.03 0.04 0.05 0.06],'xtick',[0 10 20 30 40 50 60]);
hold off

% Inset
insx=0.55;
insy=0.57;
insxl = 0.4;
insyl = 0.37;
axes('position',[insx insy insxl insyl])
hold on

R_h = 1;
N = length(Chan6a_Z);
ind = [1:6:24,25:35,36:2:51,52:6:N];
hand(1) = plot(Chan6a_Z(ind)/R_h,Chan6a_hZ(ind)/sum(Chan6a_hZ),'linestyle','none','Marker','o','MarkerFaceColor',COLOR{1},'MarkerEdgeColor',COLOR{1},'MarkerSize',ms,'linewidth',0.75);
hand(2) = plot(Wall_Z(ind)/R_h,Wall_hZ(ind)/sum(Wall_hZ),'linestyle','none','Marker','s','MarkerFaceColor','none','MarkerEdgeColor',COLOR{2},'MarkerSize',ms+1,'linewidth',0.75);
plot([5 5],[0,1],'--b','linewidth',linewidth)

xlim([0,6])
ylim([0,0.065])

XLIM = xlim;
YLIM = ylim;
xlab = '$h/R_h$';
x_dim = x_label_x_inset*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y_inset*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',small_fontsize)
ylab = '$\mathrm{P}(h)$';
x_dim = y_label_x_inset*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y_inset*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',small_fontsize,'Rotation',90)

LEGEND = {'BW','SC'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+insx*xl+0.45*insxl,dy_b+insy*yl+0.5*insyl,0.05,0.05],'fontsize',small_fontsize);
leg.ItemTokenSize = [11,1];
leg.NumColumns = 1;
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',small_fontsize-2,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',[0,1,2,3,4,5,6],'ytick',[0,0.02,0.04,0.06]);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print('rollers_BWSC_Comparison','-dpdf','-painters')

function [x_ref,y_ref] = SPLINE(x,y,fac)
N = length(x);
x_ref = linspace(x(1),x(end),fac*N);
y_ref = spline(x,y,x_ref);
end