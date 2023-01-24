clc;clear;close all
warning off
addpath('./data')
format long
sub_label_x = -0.05;
sub_label_y = 1.15 ;

y_label_x = -0.14;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5  ;
x_label_y = -0.15;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 3;
yl_in = 2.4;
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
ms = 4;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],'c',[0.47,0.67,0.19],[0.3,0.75,0.93],[0.9 0 0],'b'};

subplot('position',pos11)
hold on
%% Monopole
load effRad_monopole.mat
Ls = double(Ls); ns = double(ns); ws = double(ws);
h = Ls(1)./ns(1);
for j = 1:length(ws)
    for k = 1:length(betas)
        for l = 1:nTrials
            p = polyfit(1./Ls',M(:,j,k,l),1);
            effRh(j,k,l) = 1/h/(6*pi*eta*p(2));
        end
    end
end
Rhs = mean(effRh,3);
err = 4*std(effRh,0,3)./Rhs*100;

betaw = ws'*betas;
Rhdhw = Rhs./(h*ws');
p = polyfit(betaw(:),Rhdhw(:),10);
bax = linspace(min(betas)*ws(1),max(betas)*ws(end));
plot(bax, polyval(p,bax), '-k','linewidth',linewidth);
f_monopole = polyval(p,bax);
for j = 1:length(ws)
    hand(j) = plot(betas*ws(j),Rhs(j,:)/(h*ws(j)),'linestyle','none','Marker','o','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',ms,'linewidth',linewidth/2);
end

%% Dipole
load effRad_dipole.mat
Ls = double(Ls); ns = double(ns); ws = double(ws);
h = Ls(1)./ns(1);
for j = 1:length(ws)
    for k = 1:length(betas)
        for l = 1:nTrials
            p = polyfit(1./Ls'.^3,M(:,j,k,l),1);
            effRh(j,k,l) = 1/h/(8*pi*eta*p(2))^(1/3);
        end
    end
end
Rhs = mean(effRh,3);
err = 4*std(effRh,0,3)./Rhs*100;

betaw = ws'*betas;
Rhdhw = Rhs./(h*ws');
p = polyfit(betaw(:),Rhdhw(:),10);
bax = linspace(min(betas)*ws(1),max(betas)*ws(end));
plot(bax, polyval(p,bax), '-k','linewidth',linewidth);
f_dipole = polyval(p,bax);
for j = 1:length(ws)
    plot(betas*ws(j),Rhs(j,:)/(h*ws(j)),'linestyle','none','Marker','o','MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',ms,'linewidth',linewidth/2);
end

xlim([0,40])
ylim([0.1,0.5])
XLIM = xlim;
YLIM = ylim;
xlab = '$\beta$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$c(\beta)=R_h/(mh)$';
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

LEGEND = {'$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.6*xl,dy_b+0.75*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
leg.NumColumns = 3;
title(leg,'$m$','interpreter','latex')
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',[0.1,0.2,0.3,0.4,0.5]);
hold off

print('cbetafit','-dpdf','-painters')