clc;clear;close all
warning off
format long
sub_label_x = -0.17;
sub_label_y = 0.97 ;

y_label_x = -0.19;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.17;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 2;
dx_l_in = 0.65;dx_r_in = 0.12;
dy_b_in = 0.45;dy_t_in = 0.1;
xwidth = 2*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 2 ;
figbox_linewidth = 1;
ms = 3;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
pos12 = [(dx_l+xl)+dx_l,dy_b,xl,yl];

dat = readmatrix(['lanczos_SC_to_BW.dat']);
idx  = find(isnan(dat(:,1)) == 1);
SCBW1   = dat(1:idx(1)-1,1:2);
SCBW2   = dat(idx(1)+1:idx(2)-1,1:2);
SCBW3  = dat(idx(2)+1:idx(3)-1,1:2);
BW = dat(idx(3)+1:end,1:2);

dat = readmatrix(['lanczos_SC_to_TP.dat']);
idx  = find(isnan(dat(:,1)) == 1);
SCTP1   = dat(1:idx(1)-1,1:2);
SCTP2   = dat(idx(1)+1:idx(2)-1,1:2);
SCTP3  = dat(idx(2)+1:idx(3)-1,1:2);
SCTP4  = dat(idx(3)+1:idx(4)-1,1:2);
TP = dat(idx(4)+1:end,1:2);

% COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],[0,0.45,0.74],[0.47,0.67,0.19],[0.3,0.75,0.93]};
COLOR = {[0.9 0 0] 'b' [0,0.7,0] 'k' [0.93,0.69,0.13]};
MARKER = {'^','s','d','o','V'};
MS = ms*ones(1,length(MARKER));
MS(2) = ms+1;

%%
subplot('position',pos11)
hold on

i = 1;hand(i) = plot(SCBW1(:,1),SCBW1(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 2;hand(i) = plot(SCBW2(:,1),SCBW2(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 3;hand(i) = plot(SCBW3(:,1),SCBW3(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 4;hand(i) = plot(BW(:,1),BW(:,2),'linestyle','none','color',COLOR{i},'Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));
i = 4;plot(BW(:,1),BW(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);

xlim([0,40])
ylim([1e-8,1])
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

z0 = 3e-8;H = 25;
x0 = 2   ;L = 15;
n = 12;
x = linspace(x0,x0+L,n);
dhamp = 0.05;
plot([x0,x0+L],[z0,z0],'-k','linewidth',linewidth)
plot([x0,x0+L],H*[z0,z0],'-k','linewidth',linewidth)
loc = 0.2;
dh = (1-2*rand(1,n))*dhamp;
plot(x(2:n-1),10.^((loc+dh(2:n-1))*(log10(H*z0)-log10(z0))+log10(z0)),'ok','MarkerSize',4,'LineWidth',0.5)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'SC $H=7.5R_h$','SC $H=17.5R_h$','SC $H=39R_h$','BW'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.6*xl,dy_b+0.75*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-8,1e-6,1e-4,1e-2,1e0]);
hold off

%%
subplot('position',pos12)
hold on
COLOR = {[0.9 0 0] 'b' [0,0.7,0] [0.93,0.69,0.13] 'k'};
MARKER = {'^','s','d','V','o'};

i = 1;plot(SCTP1(:,1),SCTP1(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
hand(i) = plot(SCTP1(:,1),SCTP1(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

i = 2;plot(SCTP2(:,1),SCTP2(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
hand(i) = plot(SCTP2(:,1),SCTP2(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

i = 3;plot(SCTP3(:,1),SCTP3(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
hand(i) = plot(SCTP3(:,1),SCTP3(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

i = 4;plot(SCTP4(:,1),SCTP4(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
hand(i) = plot(SCTP4(:,1),SCTP4(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

i = 5;plot(TP(:,1),TP(:,2),'-','color',COLOR{i},'linewidth',linewidth/2);
hand(i) = plot(TP(:,1),TP(:,2),'linestyle','none','Marker',MARKER{i},'MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',MS(i));

xlim([0,40])
ylim([1e-8,1])
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

plot([x0,x0+L],[z0,z0],'-k','linewidth',linewidth)
plot([x0,x0+L],H*[z0,z0],'-k','linewidth',linewidth)
loc = 0.5;
dh = (1-2*rand(1,n))*dhamp;
plot(x(2:n-1),10.^((loc+dh(2:n-1))*(log10(H*z0)-log10(z0))+log10(z0)),'ok','MarkerSize',4,'LineWidth',0.5)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'SC $H=8.5R_h$','SC $H=17.5R_h$','SC $H=39R_h$','SC $H=107R_h$','TP ($H=400R_h$)'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[(dx_l+xl)+dx_l+0.66*xl,dy_b+0.73*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[1e-8,1e-6,1e-4,1e-2,1e0]);
hold off

print('lanczos_SC_to_BWTP','-dpdf','-painters')