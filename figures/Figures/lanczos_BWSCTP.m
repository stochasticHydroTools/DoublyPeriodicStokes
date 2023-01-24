clc;clear;close all
warning off
format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.09;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.1;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 3;
yl_in = 2.4;
dx_l_in = 0.5;dx_r_in = 0.12;
dy_b_in = 0.4;dy_t_in = 0.1;
xwidth = 1*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 2 ;
figbox_linewidth = 1;
ms = 5;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

dat = readmatrix(['lanczos_BWSCTP.dat']);
idx  = find(isnan(dat(:,1)) == 1);
BW2048   = dat(1:idx(1)-1,1:2);
BW8192   = dat(idx(1)+1:idx(2)-1,1:2);
BW32768  = dat(idx(2)+1:idx(3)-1,1:2);
BW131072 = dat(idx(3)+1:idx(4)-1,1:2);
SC2048   = dat(idx(4)+1:idx(5)-1,1:2);
SC8192   = dat(idx(5)+1:idx(6)-1,1:2);
SC32768  = dat(idx(6)+1:idx(7)-1,1:2);
SC131072 = dat(idx(7)+1:idx(8)-1,1:2);
TP2048   = dat(idx(8)+1:idx(9)-1,1:2);
TP8192   = dat(idx(9)+1:idx(10)-1,1:2);
TP32768  = dat(idx(10)+1:idx(11)-1,1:2);
TP131072 = dat(idx(11)+1:end,1:2);

%%
subplot('position',pos11)
hold on
COLOR = {[0.9 0 0] 'b' [0,0.7,0] 'k' [0.93,0.69,0.13]};
MARKER = {'^','s','d','o'};
MS = ms*ones(1,length(MARKER));
MS(2) = ms+0.75;

i = 4;plot(BW131072(:,1),BW131072(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
i = 1;hand(i) = plot(BW2048(:,1),BW2048(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 2;hand(i) = plot(BW8192(:,1),BW8192(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 3;hand(i) = plot(BW32768(:,1),BW32768(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 4;hand(i) = plot(BW131072(:,1),BW131072(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

i = 4;plot(SC131072(:,1),SC131072(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
i = 1;plot(SC2048(:,1),SC2048(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 2;plot(SC8192(:,1),SC8192(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 3;plot(SC32768(:,1),SC32768(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 4;plot(SC131072(:,1),SC131072(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));


i = 1;plot(TP2048(:,1),TP2048(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
i = 2;plot(TP8192(:,1),TP8192(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
i = 3;plot(TP32768(:,1),TP32768(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
i = 4;plot(TP131072(:,1),TP131072(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);

i = 1;plot(TP2048(:,1),TP2048(:,2),'linestyle','none','color',COLOR{i},'Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 2;plot(TP8192(:,1),TP8192(:,2),'linestyle','none','color',COLOR{i},'Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 3;plot(TP32768(:,1),TP32768(:,2),'linestyle','none','color',COLOR{i},'Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 4;plot(TP131072(:,1),TP131072(:,2),'linestyle','none','color',COLOR{i},'Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

xlim([0,30])
ylim([1e-6,1])
XLIM = xlim;
YLIM = ylim;
xlab = '$n$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\epsilon_n$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

x_dim = 22.5;
y_dim = 1e-5;
clear text
text(x_dim,y_dim,'BW','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
x_dim = 12;
y_dim = 1e-5;
clear text
text(x_dim,y_dim,'SC','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
x_dim = 20;
y_dim = 1e-2;
clear text
text(x_dim,y_dim,'TP','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'$2^{11}$','$2^{13}$','$2^{15}$','$2^{17}$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.1*xl,dy_b+0.2*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
legend boxoff
title(leg,'$N$','interpreter','latex')

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k');
hold off

print('lanczos_BWSCTP','-dpdf','-painters')