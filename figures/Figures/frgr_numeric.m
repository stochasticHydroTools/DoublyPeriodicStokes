clc;clear;close all
warning off
format long
sub_label_x = -0.1;
sub_label_y = 1   ;

y_label_x = -0.16;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.13;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 2.5;
dx_l_in = 0.57;dx_r_in = 0.1;
dy_b_in = 0.45;dy_t_in = 0.2;
xwidth = 2*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1 ;
figbox_linewidth = 1;
ms = 1;
ms_legend = 10;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
pos12 = [(dx_l+xl)+dx_l,dy_b,xl,yl];

%this script generate pair mobility function f(r) and g(r) for TP case.
%and validate through comparing with FCM analytical result (ref:maxey 2001)
addpath('./data')
transparency=0.1; %it sets the traspanrency of the scatter plots for w=6 curves

%% Calculations
% load output data from pair_perp_check_TP.py
mobx_u_w6 = load('perpendicular_mobility_unit_w6.txt');
mobx_u_w4 = load('perpendicular_mobility_unit_w4.txt');

mobx_para = load('parallel_mobility_unit_w4.txt');
mobx_u_w4_gr = mobx_para-mobx_u_w4;

mobx_para  = load('parallel_mobility_unit_w6.txt');
mobx_u_w6_gr = mobx_para-mobx_u_w6;

% set sim params for unit/non_unit grid spacing cases
eta = 1/4/sqrt(pi);
F = 1;

%1.5539 for w=6; 1.3437 for w=5; 1.2047 for w=4 
Rh_u_w6 = 1.5539;
Rh_u_w4 = 1.2047;

mu0_u_w6 = 1/(6*pi*eta*Rh_u_w6); %mobility in free space
mu0_u_w4 = 1/(6*pi*eta*Rh_u_w4); %mobility in free space

Ls_u = linspace(0,25,1000); %tried pair distance in [0, 25]

nTrials = 10; %number of trials for each distance Ls_u
Lbox = 200;

%compute the anlaytical mob function f(r) with FCM kernel as reference
sigma_w6 = Rh_u_w6*sqrt(2)/sqrt(pi); %sigma for FCM Gaussians
sigma_w4 = Rh_u_w4*sqrt(2)/sqrt(pi);

Ls_u(1) = 1e-6; %just to temporary remove singularity at zero

%generate the analytic FCM kernel data
fr_FCM=((1.0+sigma_w6*sigma_w6./Ls_u./Ls_u).*erf(Ls_u/sqrt(2.0)/sigma_w6)...
    -2*sigma_w6*(2.0*pi)^(-0.5).*exp(-Ls_u.*Ls_u/2/sigma_w6/sigma_w6)./Ls_u)...
    /8.0/pi/eta./Ls_u;

gr_FCM=((1.0-3.0*sigma_w6*sigma_w6./Ls_u./Ls_u).*erf(Ls_u/sqrt(2.0)/sigma_w6)...
    +6.0*sigma_w6*(2.0*pi)^(-0.5).*exp(-Ls_u.*Ls_u/2/sigma_w6/sigma_w6)./Ls_u)...
    /8.0/pi/eta./Ls_u;

fr_FCM_w4=((1.0+sigma_w4*sigma_w4./Ls_u./Ls_u).*erf(Ls_u/sqrt(2.0)/sigma_w4)...
    -2*sigma_w4*(2.0*pi)^(-0.5).*exp(-Ls_u.*Ls_u/2/sigma_w4/sigma_w4)./Ls_u)...
    /8.0/pi/eta./Ls_u;

gr_FCM_w4=((1.0-3.0*sigma_w4*sigma_w4./Ls_u./Ls_u).*erf(Ls_u/sqrt(2.0)/sigma_w4)...
    +6.0*sigma_w4*(2.0*pi)^(-0.5).*exp(-Ls_u.*Ls_u/2/sigma_w4/sigma_w4)./Ls_u)...
    /8.0/pi/eta./Ls_u;

% difference in pair mobility function f(r) compare with FCM ref
rel_diff_fr_w6 = (mobx_u_w6' + 2.8373*Rh_u_w6/Lbox/(6*pi*eta*Rh_u_w6) - fr_FCM);
rel_diff_fr_w4 = (mobx_u_w4' + 2.8373*Rh_u_w4/Lbox/(6*pi*eta*Rh_u_w4) - fr_FCM_w4);

rel_diff_gr_w4 = (mobx_u_w4_gr' - gr_FCM_w4);
rel_diff_gr_w6 = (mobx_u_w6_gr' - gr_FCM);

Ls_u(1) = 0; %shift back to zero

COLOR = {'c' 'r' 'k' 'g' 'm' 'k'};
%%
subplot('position',pos11)
hold on

i = 1;
plot(Ls_u'/Rh_u_w4,8*pi*eta*Ls_u'.*(mobx_u_w4+2.8373*Rh_u_w4/Lbox/(6*pi*eta*Rh_u_w4)),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 2;
plot(Ls_u'/Rh_u_w6,8*pi*eta*Ls_u'.*(mobx_u_w6+2.8373*Rh_u_w6/Lbox/(6*pi*eta*Rh_u_w6)),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 3;
plot(Ls_u/Rh_u_w6,8*pi*eta*Ls_u.*fr_FCM,'-','color',COLOR{i},'LineWidth',linewidth);
hand(i) = plot([1,100],[100,100],'-','color',COLOR{i},'LineWidth',linewidth);

i = 4;
plot(Ls_u'/Rh_u_w4,8*pi*eta*Ls_u'.*mobx_u_w4_gr,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 5;
plot(Ls_u'/Rh_u_w6,8*pi*eta*Ls_u'.*mobx_u_w6_gr,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 6;
plot(Ls_u/Rh_u_w6,8*pi*eta*Ls_u.*gr_FCM,'--','color',COLOR{i},'LineWidth',linewidth);
hand(i) = plot([1,100],[100,100],'--','color',COLOR{i},'LineWidth',linewidth);

xlim([0 15])
ylim([0 1.16])
XLIM = xlim;
YLIM = ylim;
xlab = '$d/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\left(8\pi\eta d\right)\times f(d)$ and $g(d)$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

LEGEND = {'$f(d)$, $m=4$','$f(d)$, $m=6$','$f_{\mathrm{FCM}}(d)$ reference','$g(d)$, $m=4$','$g(d)$, $m=6$','$g_{\mathrm{FCM}}(d)$ reference'};
legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.6*xl,dy_b+0.25*yl,0.05,0.05]);
% leg.ItemTokenSize = [16,1];
legend boxoff;

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = sub_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k');
hold off

%%
subplot('position',pos12)
hold on

mult = 1e3;
i = 1;
for j = 1:10
%     plot(Ls_u/Rh_u_w4,8*pi*eta*Ls_u.*rel_diff_fr_w4(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
    plot(Ls_u/Rh_u_w4,mult*8*pi*eta*Rh_u_w4.*rel_diff_fr_w4(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
end

i = 2;
for j = 1:10
%     plot(Ls_u/Rh_u_w6,8*pi*eta*Ls_u.*rel_diff_fr_w6(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
    plot(Ls_u/Rh_u_w6,mult*8*pi*eta*Rh_u_w6.*rel_diff_fr_w6(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
end

i = 4;
for j = 1:10
%     plot(Ls_u/Rh_u_w4,8*pi*eta*Ls_u.*rel_diff_gr_w4(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
    plot(Ls_u/Rh_u_w4,mult*8*pi*eta*Rh_u_w4.*rel_diff_gr_w4(j,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
end

i = 5;
for j = 1:10
%     plot(Ls_u/Rh_u_w6,8*pi*eta*Ls_u.*rel_diff_gr_w6(i,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
    plot(Ls_u/Rh_u_w6,mult*8*pi*eta*Rh_u_w6.*rel_diff_gr_w6(i,:),'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms);
end

% For legends
i = 1;
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 2;
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 4;
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

i = 5;
hand(i) = plot(100,100,'linestyle','none','Marker','.','MarkerFaceColor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',ms_legend);

hand(3) = [];

xlim([0 15])
ylim([-8 10])
XLIM = xlim;
YLIM = ylim;
xlab = '$d/R_h$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\left(8\pi\eta R_h\right)\times\Delta f(d)$ and $\Delta g(d)$';
'$\left(8\pi\eta d\right)\times f(d)$ and $g(d)$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

LEGEND = {'$\Delta f(d)$, $m=4$','$\Delta f(d)$, $m=6$','$\Delta g(d)$, $m=4$','$\Delta g(d)$, $m=6$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[(dx_l+xl)+dx_l+0.6*xl,dy_b+0.75*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
legend boxoff;

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = sub_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

exp_TXT = ['$\times 10^{-' num2str(log10(mult)) '}$'];
x_exp = 0.075;
y_exp = 1.035;
x_dim = x_exp*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_exp*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,exp_TXT,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[-8 -4 0 4 8],'yticklabels',{'$-8$','$-4$','$0$','$4$','$8$'});
hold off

print('frgr_numeric','-dpdf','-painters')
