clc;clear;close all
warning off
format long
sub_label_x = -0.15;
sub_label_y = 0.95 ;

y_label_x = -0.15;
y_label_y = 0.5 ;
y_label_x_inset = -0.3;
y_label_y_inset = 0.45;

x_label_x = 0.5 ;
x_label_y = -0.1;
x_label_x_inset = 0.5  ;
x_label_y_inset = -0.29;

xl_in = 3.8;
yl_in = 3.4;
dx_l_in = 0.87;dx_r_in = 0.1;
dy_b_in = 0.5;dy_t_in = 0.1;

dx_m_in = 0.4;
dy_bsub_in = 0.1;
ylsub_in = (yl_in-dy_bsub_in)/2;
xlsub_in = ylsub_in;

xwidth = (xl_in+dx_l_in)+dx_m_in+xlsub_in+dx_r_in;
ywidth = (yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 2 ;
figbox_linewidth = 1;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;
dx_m = dx_m_in/xwidth;
dy_bsub = dy_bsub_in/ywidth;
ylsub = ylsub_in/ywidth;
xlsub = xlsub_in/xwidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
possub1 = [dx_l+xl+dx_m,dy_b+ylsub+dy_bsub,xlsub,ylsub];
possub2 = [dx_l+xl+dx_m,dy_b,xlsub,ylsub];

%% pair mobility decay
dat = readmatrix(['pair_mobilities_decay.dat']);
idx  = isnan(dat(:,1));
dat(idx,:) = [];
N = 300;
BWParallel = dat(1:N,:);
BWPerpendicular = dat(N+1:2*N,:);
SCParallel = dat(2*N+1:3*N,:);
SCPerpendicular = dat(3*N+1:4*N,:);

dat = readmatrix(['pair_mobilities_decay_theory.dat']);
idx  = isnan(dat(:,1));
dat(idx,:) = [];

BWParallelTheory = dat(1:N,:);
BWPerpendicularTheory = dat(N+1:2*N,:);
SCTheory = dat(2*N+1:3*N,:);
%%
% COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],[0,0.45,0.74],[0.47,0.67,0.19],[0.3,0.75,0.93]};
subplot('position',pos11)
hold on
COLOR = {'k' [0.9 0 0]};
FORM = {'-',':'};

H = 8;eta = 1;R_h = 1;
d = SCParallel(:,1)*R_h;
hand(1) = plot(d/R_h,abs(BWParallel(:,2)),'linestyle',FORM{1},'color',COLOR{2},'linewidth',linewidth);
hand(2) = plot(d/R_h,abs(BWPerpendicular(:,2)),'linestyle',FORM{2},'color',COLOR{2},'linewidth',linewidth);
hand(3) = plot(d/R_h,abs(SCParallel(:,2)),'linestyle',FORM{1},'color',COLOR{1},'linewidth',linewidth);
hand(4) = plot(d/R_h,abs(SCPerpendicular(:,2)),'linestyle',FORM{2},'color',COLOR{1},'linewidth',linewidth);

% hand(5) = plot(d/R_h,BWParallelTheory(:,2),'--','color','b','linewidth',linewidth);
% hand(6) = plot(d/R_h,BWPerpendicularTheory(:,2),'--','color',[0.93,0.69,0.13],'linewidth',linewidth);
% hand(7) = plot(d/R_h,SCTheory(:,2),'--','color',[0. 0.7 0],'linewidth',linewidth);
 
mu0_BW = 0.04567698;% Calculated numerically for d=0.001R_h
mu0_SC = 0.04004688;
hand(5) = plot(d/R_h,3*H^2./(8*pi*eta*d.^3)/mu0_BW,'--','color',[0. 0.7 0],'linewidth',linewidth);
hand(6) = plot(d/R_h,3*H^4./(64*pi*eta*d.^5)/mu0_BW,'--','color',[0.5 0 0.5],'linewidth',linewidth);
hand(7) = plot(d/R_h,3*H./(32*pi*eta*d.^2)/mu0_SC,'--','color',[0.85,0.33,0.1],'linewidth',linewidth);

xlim([1e-1,128])
ylim([1e-4,2e0])
XLIM = xlim;
YLIM = ylim;
xlab = '$d/R_h$';
x_dim = 10^(x_label_x*(log10(XLIM(2))-log10(XLIM(1)))+log10(XLIM(1)));
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$\displaystyle\frac{\mid\nu^{tt}\mid}{\nu^{tt}(0)}$';
x_dim = 10^(y_label_x*(log10(XLIM(2))-log10(XLIM(1)))+log10(XLIM(1)));
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'BW $\nu_{\parallel}^{tt}$','BW $\nu_{\perp}^{tt}$','SC $\nu_{\parallel}^{tt}$','SC $\nu_{\perp}^{tt}$','$\displaystyle\frac{3H^2}{8\pi\eta d^3}$','$\displaystyle\frac{3H^4}{64\pi\eta d^5}$','$\displaystyle\frac{3H}{32\pi\eta d^2}$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.18*xl,dy_b+0.46*yl,0.05,0.3],'fontsize',small_fontsize-1);
leg.ItemTokenSize = [22,1];
leg.NumColumns = 2;
legend boxoff
% title(leg,'$\alpha$','interpreter','latex')

sub_label_y_modified = ((1+sub_label_y)*ylsub_in+dy_bsub_in)/(yl_in);
x_dim = 10^(0.7*sub_label_x*(log10(XLIM(2))-log10(XLIM(1)))+log10(XLIM(1)));
y_dim = 10^(sub_label_y_modified*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'xscale','log','yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k');%,'xtick',XTICK,'xticklabels',XTICKLABELS,'ytick',YTICK,'yticklabels',YTICKLABELS);
hold off

% Inset
insx=0.22;
insy=0.24;
axes('position',[insx insy 0.18 0.24])
hold on

hand(1) = plot(d,SCParallel(:,2),'linestyle',FORM{1},'color',COLOR{1},'linewidth',linewidth);
hand(2) = plot(d,SCPerpendicular(:,2),'linestyle',FORM{2},'color',COLOR{1},'linewidth',linewidth);
hand(3) = plot(d,BWParallel(:,2),'linestyle',FORM{1},'color',COLOR{2},'linewidth',linewidth);
hand(4) = plot(d,BWPerpendicular(:,2),'linestyle',FORM{2},'color',COLOR{2},'linewidth',linewidth);

xlim([0,25])
ylim([-0.07,0.1])
XLIM = xlim;
YLIM = ylim;
xlab = '$d/R_h$';
x_dim = x_label_x_inset*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y_inset*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)
ylab = '$\displaystyle\frac{\nu^{tt}}{\nu^{tt}(0)}$';
x_dim = y_label_x_inset*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y_inset*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)

box on
set(gca,'linewidth',figbox_linewidth/2,'fontsize',small_fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',[0 10 20],'xticklabel',{'$0$','$10$','$20$'},'ytick',[-0.05 0 0.05 0.1],'yticklabel',{'$-0.05$','$0$','$0.05$','$0.1$'});
hold off

%% flow pattern
subplotpos = [possub1;possub2];
DATA = {'flowpattern_BW.dat','flowpattern_SC.dat'};
sublabel = {'(b)','(c)'};
geometry = {'BW','SC'};
for example = 1:2
    subplot('position',subplotpos(example,:))
    hold on

dat = readmatrix(DATA{example});

[x,y,u,v,N] = datExtract(dat);
[psi] = streamlineCalc(u,v,x,y,N,1);
psi_Mirr = psi(:,end:-1:1);
[xx,yy] = meshgrid(x,y);

refine_ratio = 4;% must be even
x_fine = linspace(x(1),x(end),refine_ratio*length(x));
y_fine = linspace(y(1),y(end),refine_ratio*length(y));
[xx_fine,yy_fine] = meshgrid(x_fine,y_fine);
psi_fine = interp2(xx,yy,psi,xx_fine,yy_fine,'cubic');
psi_fine_Mirr = psi_fine(:,end:-1:1);

levels = 6;
% contour(xx(N/2+1:end,:),yy(N/2+1:end,:),psi(:,N/2+1:end)',levels,'color','k');shading flat;
% contour(xx(1:N/2,:),yy(1:N/2,:),psi_Mirr(:,1:N/2)',levels,'color','k');shading flat;
% contour(xx,yy,psi',2*levels+1,'color','r','linestyle','--');shading flat;colormap(cm);
N_fine = length(x_fine);
contour(xx_fine(N_fine/2+1:end,:),yy_fine(N_fine/2+1:end,:),psi_fine(:,N_fine/2+1:end)',levels,'color','k');shading flat;
contour(xx_fine(1:N_fine/2,:),yy_fine(1:N_fine/2,:),psi_fine_Mirr(:,1:N_fine/2)',levels,'color','k');shading flat;


threshold = 1;
u = u.*(abs(xx.^2+yy.^2) > threshold);
v = v.*(abs(xx.^2+yy.^2) > threshold);

%% Force arrow
HeadWidth  = 7 ;
HeadLength = 7 ;
LineLength = 4 ;
LineWidth = 1  ;
COLOR = 'r';
ah_r = annotation('arrow','headStyle','plain','HeadLength',HeadLength,'HeadWidth',HeadWidth,'color',COLOR,'linewidth',LineWidth);
set(ah_r,'parent',gca);
set(ah_r,'position',[1 0 LineLength 0]);

F_X = 1+LineLength+1.3;F_Y = 0;
clear text
text(F_X,F_Y,'$F_x$','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

%% Velocity Arrows
%% at select points
points = [5,2;-5,2];[pM,~] = size(points);
vx = zeros(1,pM);vy = zeros(1,pM);
vel_value = zeros(1,pM);
for i = 1:pM
    i_x = find(x >= points(i,1),1);
    i_y = find(y >= points(i,2),1);
    vx(i) = u(i_x,i_y);
    vy(i) = v(i_x,i_y);
    vel_value(i) = sqrt(vx(i)^2+vy(i)^2);
end

HeadWidth  = 4  ;
HeadLength = 4  ;
LineLength = 1.5;
LineWidth = 0.75;
for i = 1:pM
    COLOR = 'b';
    ah_r = annotation('arrow','headStyle','plain','HeadLength',HeadLength,'HeadWidth',HeadWidth,'color',COLOR,'linewidth',LineWidth);
    set(ah_r,'parent',gca);
    set(ah_r,'position',[points(i,1) points(i,2) LineLength*vx(i)/vel_value(i) LineLength*vy(i)/vel_value(i)]);

    if points(i,2) > 0.1
        ah_r = annotation('arrow','headStyle','plain','HeadLength',HeadLength,'HeadWidth',HeadWidth,'color',COLOR,'linewidth',LineWidth);
        set(ah_r,'parent',gca);
        set(ah_r,'position',[points(i,1) -points(i,2) LineLength*vx(i)/vel_value(i) -LineLength*vy(i)/vel_value(i)]);
    end

end

%% sphere
rad = 1;
cent = [0 0];
pos = [cent-rad 2*rad 2*rad];
rect = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0.72 0.45 0.2], 'Edgecolor','none');

xlim([-11,11])
ylim([-11,11])
XLIM = xlim;
YLIM = ylim;
if example == 2
    xlab = '$x/R_h$';
    x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = (yl/ylsub)*x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
end
ylab = '$y/R_h$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

INFO_X = 0.05;INFO_Y = 0.5;
x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = INFO_Y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,geometry{example},'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','latex','fontsize',fontsize)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = sub_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,sublabel{example},'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

ytick = -8:4:8;
yticklabels = cell(1,length(ytick));
for i = 1:length(ytick)
    yticklabels{i} = ['$' num2str(ytick(i)) '$'];
end
xtick = [];
if example == 2
    xtick = ytick;
    xticklabels = yticklabels;
end

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',xtick,'xticklabels',xticklabels,'ytick',ytick,'yticklabels',yticklabels);
hold off

end

print('pair_mobilities_decay_flowpattern','-dpdf','-painters')

function [x,y,u,v,N] = datExtract(dat)
N = sqrt(length(dat));
y = dat(1:N,2)';x = y;
u = zeros(N,N);v = zeros(N,N);
for i = 1:N
    for j = 1:N
        u(i,j) = dat((i-1)*N+j,3);
        v(i,j) = dat((i-1)*N+j,4);
    end
end
end

function [psi] = streamlineCalc(u,v,x,y,N,parallel)
psi = zeros(N,N);
if parallel
    u_middle = 1/2*(u(:,N/2)+u(:,N/2+1));
    for i = 1:N
        for j = 1:N
            if j > N/2
                y_integration = [0,y(N/2+1:j)];
                u_integration = [u_middle(i),u(i,N/2+1:j)];
            else
                y_integration = [0,y(N/2:-1:j)];
                u_integration = [u_middle(i),u(i,N/2:-1:j)];
            end
            psi(i,j) = trapz(y_integration,u_integration);
        end
    end
else
    v_middle = 1/2*(v(N/2,:)+v(N/2+1,:));
    for i = 1:N
        for j = 1:N
            if i > N/2
                x_integration = [0,x(N/2+1:i)];
                v_integration = [v_middle(i),v(N/2+1:i,j)'];
            else
                x_integration = [0,x(N/2:-1:i)];
                v_integration = [v_middle(i),v(N/2:-1:i,j)'];
            end
            psi(i,j) = trapz(x_integration,v_integration);
        end
    end
end

end