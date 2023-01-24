clc;clear;close all
warning off
addpath('./data')

format long

sub_label_x = -0.22;
sub_label_y = 0.97 ;

y_label_x = -0.22;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.17;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 2;
dx_l_in = 0.75;dx_r_in = 0.12;
dy_b_in = 0.45;dy_t_in = 0.1;
xwidth = 2*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1.5 ;
figbox_linewidth = 1;
ms = 5;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
pos12 = [(dx_l+xl)+dx_l,dy_b,xl,yl];

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
%% Monopole
subplot('position',pos11)
hold on

LEGEND = cell(1,4);
hand(4) = plot(H_bw_ref_noimg, abs(M_bw_ref_noimg(:,2,4)),'--','color',[0.9 0 0],'linewidth',linewidth);
LEGEND{4} = 'RPB';
hand(3) = plot(H_bw_ref, abs(M_bw_ref(:,2,4)),'-','color',[0 0.7 0],'linewidth',linewidth);
LEGEND{3} = 'PRPB';
hand(1) = plot(H_bw, abs(M_bw(:,2,4)),'s','MarkerFaceColor',COLOR{3},'MarkerEdgeColor',COLOR{3},'MarkerSize',ms+1,'linewidth',0.5);
LEGEND{1} = '$m\!=\!6$';
hand(2) = plot(H_bw_ref2x, abs(M_bw_ref2x(:,2,4)),'o-','color',COLOR{4},'MarkerFaceColor',COLOR{4},'MarkerSize', ms-2,'linewidth',0.5);
LEGEND{2} = '$m\!=\!12$';

xlim([0,9.6])
ylim([1e-5,1e-1])
XLIM = xlim;
YLIM = ylim;
xlab = '$z/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\tilde{\mu}_{yx}^{tr}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

INFO_X = 0.3;
INFO_Y = 0.9;
x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'BW','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

legx = dx_l+0.7*xl;
legy = dy_b+0.775*yl;
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[legx,legy,0.05,0.05]);
leg.NumColumns = 1;
% leg.ItemTokenSize = [20,1];
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-5,1e-4,1e-3,1e-2,1e-1]);
hold off

%% Dipole
subplot('position',pos12)
hold on

LEGEND = cell(1,2);
hand(1) = plot(H_sc, abs(M_sc(:,2,4)),'s','MarkerFaceColor',COLOR{3},'MarkerEdgeColor',COLOR{3},'MarkerSize',ms+1,'linewidth',0.5);
LEGEND{1} = '$m\!=\!6$';
hand(2) = plot(H_sc_ref, abs(M_sc_ref(:,2,4)),'o-','color',COLOR{4},'MarkerFaceColor',COLOR{4},'MarkerSize', ms-2,'linewidth',0.5);
LEGEND{2} = '$m\!=\!12$';

xlim([0,9.6])
ylim([1e-5,1e-1])
XLIM = xlim;
YLIM = ylim;
xlab = '$z/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\mid\!\!\tilde{\mu}_{yx}^{tr}\!\!\mid$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)


x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'SC','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

legx = xl+2*dx_l+0.7*xl;
legy = dy_b+0.85*yl;
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[legx,legy,0.05,0.05]);
leg.NumColumns = 1;
% leg.ItemTokenSize = [20,1];
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-5,1e-4,1e-3,1e-2,1e-1]);
hold off

print('self_tr','-dpdf','-painters')