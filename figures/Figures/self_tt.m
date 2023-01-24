clc;clear;close all
warning off
addpath('./data')
format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.10;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.13;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 3;
yl_in = 2.4;
dx_l_in = 0.53;dx_r_in = 0.1;
dy_b_in = 0.5;dy_t_in = 0.1;
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
MS = ms*ones(1,10);MS([3,4]) = ms+1;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

load self_mobility_bw.mat; H_bw = heights; M_bw = M; 
load self_mobility_bw_ref.mat; H_bw_ref = heights; M_bw_ref = M;
load self_mobility_bw_ref_noimg.mat; H_bw_ref_noimg = heights; M_bw_ref_noimg = M;
load self_mobility_bw_w4.mat; H_bw_w4 = heights; M_bw_w4 = M;
load self_mobility_bw_ref2x.mat; H_bw_ref2x = heights; M_bw_ref2x = abs(M);

load self_mobility_sc.mat; H_sc = heights; M_sc = M;
load self_mobility_sc_ref.mat; H_sc_ref = heights; M_sc_ref = M;
load self_mobility_sc_w4.mat; H_sc_w4 = heights; M_sc_w4 = M;
data = dlmread('mobilitySlitChannel.width.19.2Rh.1blob.dat');

COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.93,0.69,0.13],[0.49,0.18,0.56],'c',[0.47,0.67,0.19],[0.3,0.75,0.93],[0.9 0 0],'b'};

subplot('position',pos11)
hold on
%% trans-trans
L = 19.2;
H = L/2;
faxen_half = (1-1.004/H+0.418*(1/H)^3+0.21*(1/H)^4-0.169*(1/H)^5);
H = L/4;
faxen_quarter = (1-0.6526/H+0.1475*(1/H)^3-0.131*(1/H)^4-0.0644*(1/H)^5);

LEGEND = cell(1,13);

% BW
j = 3;
hand(10) = plot(H_bw_ref, M_bw_ref(:,1,1),'-','color',COLOR{j},'linewidth',linewidth);
LEGEND{10} = '$\mu_{xx}^{tt}$ PRPB';
hand(8) = plot(H_bw, M_bw(:,1,1),'s','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{8} = '$\mu_{xx}^{tt}$ $m=6$';
hand(9) = plot(H_bw_w4, M_bw_w4(:,1,1),'s','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
LEGEND{9} = '$\mu_{xx}^{tt}$ $m=4$';
j = 4;
hand(13) = plot(H_bw_ref, M_bw_ref(:,3,3),'-','color', COLOR{j},'linewidth',linewidth);
LEGEND{13} = '$\mu_{zz}^{tt}$ PRPB';
hand(11) = plot(H_bw, M_bw(:,3,3),'s','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{11} = '$\mu_{zz}^{tt}$ $m=6$';
hand(12) = plot(H_bw_w4, M_bw_w4(:,3,3),'s','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
LEGEND{12} = '$\mu_{zz}^{tt}$ $m=4$';

% SC
j = 1;
hand(3) = plot(data(1:44,1),data(1:44,2),'-','color',COLOR{j},'linewidth',linewidth);
LEGEND{3} = '$\mu_{xx}^{tt}$ IBM';
hand(1) = plot(H_sc, M_sc(:,1,1),'o','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$\mu_{xx}^{tt}$ $m=6$';
hand(2) = plot(H_sc_w4, M_sc_w4(:,1,1),'o','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
LEGEND{2} = '$\mu_{xx}^{tt}$ $m=4$';
hand(7) = plot(L/4,faxen_quarter,'ko', 'MarkerFaceColor','c','MarkerSize',ms,'linewidth',0.5);
plot(L/2,faxen_half,'ko','MarkerFaceColor','c','MarkerSize',ms,'linewidth',0.5);
LEGEND{7} = '$\mu_{xx}^{tt}$ Faxen';
j = 2;
hand(6) = plot(data(1:44,1),data(46:89,2),'-','color',COLOR{j},'linewidth',linewidth);
LEGEND{6} = '$\mu_{zz}^{tt}$ IBM';
hand(4) = plot(H_sc, M_sc(:,3,3),'o','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{4} = '$\mu_{zz}^{tt}$ $m=6$';
hand(5) = plot(H_sc_w4, M_sc_w4(:,3,3),'o','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j)-2,'linewidth',0.5);
LEGEND{5} = '$\mu_{zz}^{tt}$ $m=4$';


xlim([0,9.6])
ylim([0,1])
XLIM = xlim;
YLIM = ylim;
xlab = '$z/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\tilde{\mu}^{tt}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

legx = dx_l+0.6*xl;
legy = dy_b+0.27*yl;
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[legx,legy,0.05,0.05]);
leg.ItemTokenSize = [16,1];
leg.NumColumns = 2;
legend boxoff

legdx = 0.25;
legdy = 0.45;
clear text
text(legx-0.5*legdx,legy+legdy/2,'SC','units','normalized','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)
clear text
text(legx+0.79*legdx,legy+legdy/2,'BW','units','normalized','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[0,0.2,0.4,0.6,0.8,1]);
hold off

print('self_tt','-dpdf','-painters')