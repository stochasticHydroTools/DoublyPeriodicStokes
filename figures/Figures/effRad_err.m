clc;clear;close all
warning off
addpath('./data')
format long
sub_label_x = 0.05;
sub_label_y = 0.95;

y_label_x = -0.19;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.14;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 2.25;
dx_l_in = 0.65;dx_r_in = 0.15;dx_m_in = 0.15;
dy_b_in = 0.45;dy_t_in = 0.1;
xwidth = 2*xl_in+dx_l_in+dx_r_in+dx_m_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1.5 ;
figbox_linewidth = 1;
ms = 5;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;dx_m = dx_m_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
pos12 = [(dx_l+xl)+dx_m,dy_b,xl,yl];

COLOR = {[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],'c',[0.47,0.67,0.19],[0.3,0.75,0.93],[0.9 0 0],'b'};
MARKER = {'x','s','^','V','+','d','>','<','o','p','h'};
MS = ms*ones(1,10);MS([1,2,5]) = ms+1;

%% Monopole
subplot('position',pos11)
hold on
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

bax = linspace(min(betas),max(betas));
for j = 1:length(ws)
    splnE(j) = pchip(betas*ws(j),err(j,:));
    hand(j) = plot(betas,err(j,:),'linestyle','none','Marker',MARKER{j},'MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',linewidth/2);
    bax_err = linspace(min(betas*ws(j)),max(betas)*ws(j),1000);
    splneval = ppval(splnE(j),bax_err);
    plot(bax,ppval(pchip(betas,err(j,:)),bax),'-','color',COLOR{j},'linewidth',linewidth/2);
end

xlim([1,3])
ylim([1e-8,1e2])
XLIM = xlim;
YLIM = ylim;
xlab = '$\beta/m$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '\%-error';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

INFO_X = 0.5;
INFO_Y = 0.9;
x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'Monopole','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

LEGEND = {'$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.7*xl,dy_b+0.15*yl,0.05,0.05]);
leg.ItemTokenSize = [13,1];
leg.NumColumns = 3;
title(leg,'$m$','interpreter','latex')
legend boxoff
YTICK = [1e-8,1e-6,1e-4,1e-2,1e0,1e2];
box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',YTICK);
hold off

%% Dipole
subplot('position',pos12)
hold on
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

bax = linspace(min(betas),max(betas));
for j = 1:length(ws)
    splnE(j) = pchip(betas*ws(j),err(j,:));
    hand(j) = plot(betas,err(j,:),'linestyle','none','Marker',MARKER{j},'MarkerFaceColor','none','MarkerEdgeColor',COLOR{j},'MarkerSize',MS(j),'linewidth',linewidth/2);
    bax_err = linspace(min(betas*ws(j)),max(betas)*ws(j),1000);
    splneval = ppval(splnE(j),bax_err);
    plot(bax,ppval(pchip(betas,err(j,:)),bax),'-','color',COLOR{j},'linewidth',linewidth/2);
end

xlim([1,3])
ylim([1e-8,1e2])
XLIM = xlim;
YLIM = ylim;
xlab = '$\beta/m$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(x_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
% ylab = '\%-error';
% x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
% y_dim = 10^(y_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
% clear text
% text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

INFO_X = 0.5;
INFO_Y = 0.9;
x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(INFO_Y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'Dipole','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = 10^(sub_label_y*(log10(YLIM(2))-log10(YLIM(1)))+log10(YLIM(1)));
clear text
text(x_dim,y_dim,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

% LEGEND = {'$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$'};
% leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.7*xl,dy_b+0.15*yl,0.05,0.05]);
% leg.ItemTokenSize = [13,1];
% leg.NumColumns = 3;
% title(leg,'$m$','interpreter','latex')
% legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'yscale','log','TickLabelInterpreter','latex','xcolor','k','ycolor','k','ytick',YTICK,'yticklabels',{});
hold off


print('effRad_err','-dpdf','-painters')