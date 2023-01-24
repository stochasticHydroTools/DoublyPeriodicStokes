clc;clear;close all
warning off
addpath('./data')

format long
sub_label_x = 0.2;
sub_label_y = 0.95;

y_label_x = -0.165;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.3;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 0.8;
dx_l_in = 0.53;dx_r_in = 0.1;dx_m_in = 0.15;
dy_b_in = 0.35;dy_t_in = 0.1;dy_m_in = 0.1 ;
xwidth = 2*xl_in+dx_l_in+dx_r_in+dx_m_in;
ywidth = 4*yl_in+dy_b_in+dy_t_in+3*dy_m_in;

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 11;
fontsize     = 9;
small_fontsize = 8;
linewidth    = 1 ;
figbox_linewidth = 0.5;
ms = 5;
MS = [ms+1 ms ms];
MARKER = {'s','^','o'};

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;dx_m = dx_m_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;dy_m = dy_m_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b+3*(yl+dy_m),xl,yl];
pos21 = [dx_l,dy_b+2*(yl+dy_m),xl,yl];
pos31 = [dx_l,dy_b+1*(yl+dy_m),xl,yl];
pos41 = [dx_l,dy_b+0*(yl+dy_m),xl,yl];
pos12 = [dx_l+xl+dx_m,dy_b+3*(yl+dy_m),xl,yl];
pos22 = [dx_l+xl+dx_m,dy_b+2*(yl+dy_m),xl,yl];
pos32 = [dx_l+xl+dx_m,dy_b+1*(yl+dy_m),xl,yl];
pos42 = [dx_l+xl+dx_m,dy_b+0*(yl+dy_m),xl,yl];

load pair_mobility_bw.mat; H_bw = heights; M_bw = abs(M);
load pair_mobility_bw_w4.mat; H_bw_w4 = heights; M_bw_w4 = abs(M);
load pair_mobility_bw_ref.mat; H_bw_ref = heights; M_bw_ref = abs(M);
load pair_mobility_bw_ref_noimg.mat; H_bw_ref_noimg = heights; M_bw_ref_noimg = abs(M);
load mobility_asymm_posdef_bw_w4.mat; Masym_bw_w4 = Masym; Posdef_bw_w4 = Posdef;
load mobility_asymm_posdef_bw_w6.mat; Masym_bw_w6 = Masym; Posdef_bw_w6 = Posdef;

load pair_mobility_sc.mat; H_sc = heights; M_sc = abs(M);
load pair_mobility_sc_w4.mat; H_sc_w4 = heights; M_sc_w4 = abs(M);
load pair_mobility_sc_ref.mat; H_sc_ref = heights; M_sc_ref = abs(M);
load mobility_asymm_posdef_sc_w4.mat; Masym_sc_w4 = Masym; Posdef_sc_w4 = Posdef;
load mobility_asymm_posdef_sc_w6.mat; Masym_sc_w6 = Masym; Posdef_sc_w6 = Posdef;

COLOR = {'b' [0.9 0 0] [0 0.7 0]};

%% xx BW
subplot('position',pos11)
hold on
indx = 10;
indy = 7;

j = 1;plot(H_bw_ref, M_bw_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_bw_ref, M_bw_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_bw_ref, M_bw_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_bw, M_bw(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_bw_w4, M_bw_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_bw, M_bw(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_bw_w4, M_bw_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_bw, M_bw(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_bw_w4, M_bw_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim([0,9.6])
ylim([1e-6,0.05])
XLIM = xlim;
YLIM = ylim;
ylab = '$\tilde{\nu}_{xx}^{rr}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

INFO_X = 0.5;
INFO_Y = 0.15;
x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'BW','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-6,1e-5,1e-4,1e-3,1e-2],'xtick',[]);
hold off

%% xx SC
subplot('position',pos12)
hold on

j = 1;plot(H_sc_ref, M_sc_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_sc_ref, M_sc_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_sc_ref, M_sc_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_sc, M_sc(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_sc_w4, M_sc_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_sc, M_sc(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_sc_w4, M_sc_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_sc, M_sc(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_sc_w4, M_sc_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim(XLIM)
ylim(YLIM)

x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'SC','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[],'xtick',[]);
hold off

%% yy BW
subplot('position',pos21)
hold on
indx = 11;
indy = 8;

j = 1;plot(H_bw_ref, M_bw_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_bw_ref, M_bw_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_bw_ref, M_bw_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_bw, M_bw(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_bw_w4, M_bw_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_bw, M_bw(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_bw_w4, M_bw_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_bw, M_bw(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_bw_w4, M_bw_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim([0,9.6])
ylim([2e-5,0.07])
XLIM = xlim;
YLIM = ylim;
ylab = '$\tilde{\nu}_{yy}^{rr}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-5,1e-4,1e-3,1e-2,1e-1],'xtick',[]);
hold off

%% yy SC
subplot('position',pos22)
hold on

j = 1;plot(H_sc_ref, M_sc_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_sc_ref, M_sc_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_sc_ref, M_sc_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_sc, M_sc(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_sc_w4, M_sc_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_sc, M_sc(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_sc_w4, M_sc_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_sc, M_sc(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_sc_w4, M_sc_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim(XLIM)
ylim(YLIM)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[],'xtick',[]);
hold off

%% zz BW
subplot('position',pos31)
hold on
indx = 12;
indy = 9;

j = 1;plot(H_bw_ref, M_bw_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_bw_ref, M_bw_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_bw_ref, M_bw_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_bw, M_bw(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_bw_w4, M_bw_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_bw, M_bw(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_bw_w4, M_bw_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_bw, M_bw(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_bw_w4, M_bw_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim([0,9.6])
ylim([1e-6,0.05])
XLIM = xlim;
YLIM = ylim;
ylab = '$\tilde{\nu}_{zz}^{rr}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-6,1e-5,1e-4,1e-3,1e-2],'xtick',[]);
hold off

%% zz SC
subplot('position',pos32)
hold on

j = 1;plot(H_sc_ref, M_sc_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_sc_ref, M_sc_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_sc_ref, M_sc_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_sc, M_sc(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_sc_w4, M_sc_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_sc, M_sc(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_sc_w4, M_sc_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_sc, M_sc(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_sc_w4, M_sc_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim(XLIM)
ylim(YLIM)

legx_normal = 0.5 ;
legy_normal = 0.15;
legx = dx_l+1*(xl+dx_m)+legx_normal*xl;
legy = dy_b+1*(yl+dy_m)+legy_normal*yl;
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[legx,legy,0.05,0.05]);
leg.ItemTokenSize = [4,1];
leg.NumColumns = 3;
title(leg,'$m,\;d$')
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[],'xtick',[]);
hold off

%% xz BW
subplot('position',pos41)
hold on
indx = 10;
indy = 9;

j = 1;plot(H_bw_ref, M_bw_ref(1,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_bw_ref, M_bw_ref(2,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_bw_ref, M_bw_ref(3,:,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_bw, M_bw(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_bw_w4, M_bw_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_bw, M_bw(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_bw_w4, M_bw_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_bw, M_bw(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_bw_w4, M_bw_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim([0,9.6])
ylim([2e-6,0.05])
XLIM = xlim;
YLIM = ylim;
xlim(XLIM)
ylim(YLIM)
xlab = '$z/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\tilde{\nu}_{xz}^{rr}$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-6,1e-5,1e-4,1e-3,1e-2],'xtick',[0 2 4 6 8]);
hold off

%% xz SC
subplot('position',pos42)
hold on

j = 1;plot(H_sc_ref(H_sc_ref<=9), M_sc_ref(1,H_sc_ref<=9,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 2;plot(H_sc_ref(H_sc_ref<=9), M_sc_ref(2,H_sc_ref<=9,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);
j = 3;plot(H_sc_ref(H_sc_ref<=9), M_sc_ref(3,H_sc_ref<=9,indx,indy),'-','color',COLOR{j},'linewidth',linewidth);

j = 1;
hand(1) = plot(H_sc(H_sc<=9), M_sc(1,H_sc<=9,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{1} = '$6,\;3R_h$';
% hand(4) = plot(H_sc_w4, M_sc_w4(1,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{4} = '$4,\;3R_h$';

j = 2;
hand(2) = plot(H_sc(H_sc<=9), M_sc(2,H_sc<=9,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{2} = '$6,\;4R_h$';
% hand(5) = plot(H_sc_w4, M_sc_w4(2,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{5} = '$4,\;4R_h$';

j = 3;
hand(3) = plot(H_sc(H_sc<=9), M_sc(3,H_sc<=9,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',0.5);
LEGEND{3} = '$6,\;8R_h$';
% hand(6) = plot(H_sc_w4, M_sc_w4(3,:,indx,indy),'Marker',MARKER{j},'linestyle','none','MarkerFaceColor',COLOR{j},'MarkerEdgeColor',COLOR{j}, 'MarkerSize',MS(j)-2,'linewidth',0.5);
% LEGEND{6} = '$4,\;8R_h$';

xlim(XLIM)
ylim(YLIM)
xlab = '$z/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[],'xtick',[0 2 4 6 8]);
hold off


print('pair_rr','-dpdf','-painters')