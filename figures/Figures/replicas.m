clc;clear;close all
warning off
format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.12;
y_label_y = 0.5 ;
y_label_x_inset = -0.3;
y_label_y_inset = 0.5 ;

x_label_x = 0.5 ;
x_label_y = -0.1;
x_label_x_inset = 0.5  ;
x_label_y_inset = -0.33;

xl_in = 3;
yl_in = 2.4;
dx_l_in = 0.47;dx_r_in = 0.1;
dy_b_in = 0.4;dy_t_in = 0.1;
xwidth = 1*(xl_in+dx_l_in+dx_r_in);
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

dat = readmatrix(['replicas.dat']);
idx  = find(isnan(dat(:,1)) == 1);
TitanVH9       = dat(1:idx(1)-1,1:2);
V100H16a       = dat(idx(1)+1:idx(2)-1,1:2);
V100H9aTorques = dat(idx(2)+1:idx(3)-1,1:2);
V100H9a        = dat(idx(3)+1:idx(4)-1,1:2);
FIT            = dat(idx(4)+1:end,1:2);

%%
% COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],[0,0.45,0.74],[0.47,0.67,0.19],[0.3,0.75,0.93]};
subplot('position',pos11)
hold on

COLOR = {[0,0.7,0] [0.9 0 0] 'b' [0.93,0.69,0.13] 'k'};
MARKER = {'s','o','^'};
MS = ms*ones(1,length(MARKER));
MS(1) = ms+0.75;

hand(4) = plot(FIT(:,1),FIT(:,2),'-k','linewidth',linewidth);
i = 2;hand(i) = plot(TitanVH9(:,1),TitanVH9(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
% i = 2;hand(i) = plot(V100H16a(:,1),V100H16a(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 1;hand(i) = plot(V100H9aTorques(:,1),V100H9aTorques(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 3;hand(i) = plot(V100H9a(:,1),V100H9a(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

xlim([1e3,1e6])
ylim([1,500])
XLIM = xlim;
YLIM = ylim;
xlab = '$N$';
x_dim = 10^(x_label_x*(log10(XLIM(2))-log10(XLIM(1)))+log10(XLIM(1)));
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$t\;(\mathrm{ms})$';
x_dim = 10^(y_label_x*(log10(XLIM(2))-log10(XLIM(1)))+log10(XLIM(1)));
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

LEGEND = {'V100 ($M\;\&\;D$)','Titan V ($M$)','V100 ($M$)','$\displaystyle t\!=\!5.9\!\times\!10^{-4}N$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.20*xl,dy_b+0.72*yl,0.05,0.15]);
leg.ItemTokenSize = [20,1];
legend boxoff
% title(leg,'$\alpha$','interpreter','latex')

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'xscale','log','yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k');%,'xtick',XTICK,'xticklabels',XTICKLABELS,'ytick',YTICK,'yticklabels',YTICKLABELS);
hold off

% Inset
insx=0.68;
insy=0.26;
axes('position',[insx insy 0.27 0.25])
hold on

dat = readmatrix(['replicas_initdensity.dat']);
plot(dat(:,1),dat(:,2)/sum(dat(:,2)),'linestyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);

xlim([0,6.5])
ylim([0,0.3])
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
box on
set(gca,'linewidth',figbox_linewidth/2,'fontsize',small_fontsize-2,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',[0,1,2,3,4,5,6],'ytick',[0,0.1,0.2,0.3]);
hold off

print('replicas','-dpdf','-painters')