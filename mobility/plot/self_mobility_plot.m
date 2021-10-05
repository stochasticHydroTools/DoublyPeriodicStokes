%% plot self mobility in either bottom wall or slit channel
clear all; close all; clc;
addpath('../data')
set(groot, 'defaultLineLineWidth', 2.2);
set(groot,'defaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',30);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',30);

load self_mobility_bw.mat; H_bw = heights; M_bw = M; 
load self_mobility_bw_ref.mat; H_bw_ref = heights; M_bw_ref = M;
load self_mobility_bw_ref_noimg.mat; H_bw_ref_noimg = heights; M_bw_ref_noimg = M;
load self_mobility_bw_w4.mat; H_bw_w4 = heights; M_bw_w4 = M;
load self_mobility_bw_ref2x.mat; H_bw_ref2x = heights; M_bw_ref2x = abs(M);

load self_mobility_sc.mat; H_sc = heights; M_sc = M;
load self_mobility_sc_ref.mat; H_sc_ref = heights; M_sc_ref = M;
load self_mobility_sc_w4.mat; H_sc_w4 = heights; M_sc_w4 = M;
data = dlmread('mobilitySlitChannel.width.19.2Rh.1blob.dat');
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.02], [0.1 0.05], [0.1 0.05]);


cols =  distinguishable_colors(20);

%% trans-trans
figure(1);
subplot(1,1,1);
plot(H_bw, M_bw(:,1,1),'s','color',cols(1,:), 'displayname', '$\mu_{xx}^{tt}$ $w=6$, BW'); hold on;
plot(H_bw_w4, M_bw_w4(:,1,1), 's', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname', '$\mu_{xx}^{tt}$ $w=4$, BW');
plot(H_bw_ref, M_bw_ref(:,1,1),'-','color',cols(1,:),'displayname', '$\mu_{xx}^{tt}$ PRPB, BW');
plot(H_bw, M_bw(:,3,3),'s', 'color', cols(2,:), 'displayname', '$\mu_{zz}^{tt}$ $w=6$, BW'); 
plot(H_bw_w4, M_bw_w4(:,3,3), 's', 'color', cols(2,:), 'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname', '$\mu_{zz}^{tt}$ $w=4$, BW');
plot(H_bw_ref, M_bw_ref(:,3,3),'-','color', cols(2,:),'displayname', '$\mu_{zz}^{tt}$ PRPB, BW');

L = 19.2;
H = L/2;
faxen_half = (1-1.004/H+0.418*(1/H)^3+0.21*(1/H)^4-0.169*(1/H)^5);
H = L/4;
faxen_quarter = (1-0.6526/H+0.1475*(1/H)^3-0.131*(1/H)^4-0.0644*(1/H)^5);

plot(H_sc, M_sc(:,1,1),'o','color',cols(6,:), 'displayname', '$\mu_{xx}^{tt}$ $w=6$, SC','MarkerSize', 11); hold on;
plot(H_sc_w4, M_sc_w4(:,1,1), 'o', 'color', cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname', '$\mu_{xx}^{tt}$ $w=4$, SC');
plot(data(1:44,1),data(1:44,2),'-','color', cols(6,:), 'Displayname','$\mu_{xx}^{tt}$ IBM, SC');
plot(L/2,faxen_half,'ko','MarkerFaceColor','c','MarkerSize',10,'displayname','$\mu_{xx}^{tt}$ Faxen');
plot(H_sc, M_sc(:,3,3),'o','color', cols(7,:),'displayname', '$\mu_{zz}^{tt}$ $w=6$, SC','MarkerSize', 11); 
plot(H_sc_w4, M_sc_w4(:,3,3), 'o', 'color', cols(7,:),'MarkerFaceColor',cols(7,:), 'MarkerSize',6,'displayname', '$\mu_{zz}^{tt}$ $w=4$, SC','MarkerSize', 7);
pl = plot(L/4,faxen_quarter,'ko', 'MarkerFaceColor','c','MarkerSize',10);
plot(data(1:44,1),data(46:89,2),'-','color', cols(7,:)	,'displayname','$\mu_{zz}^{tt}$ IBM, SC');
%plot(H_sc_ref, M_sc_ref(:,3,3),'-','color', cols(7,:),'displayname', '$\mu_{zz}^{tt}$ $w=12$, SC'); 

set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlabel('$z/R_h$');
ylabel('$(6\pi\eta R_h) \mu^{tt}$')
legend show;
legend boxoff;

%% rot-rot
figure(2);
subplot(1,1,1);
plot(H_bw, M_bw(:,4,4),'s','color',cols(1,:),'displayname', '$\mu_{xx}^{rr}$, BW'); hold on;
plot(H_bw_ref, M_bw_ref(:,4,4),'-','color',cols(1,:),'displayname', '$\mu_{xx}^{rr}$ PRPB, BW','linewidth',5);
plot(H_bw, M_bw(:,6,6),'s','color',cols(2,:), 'displayname', '$\mu_{zz}^{rr}$, BW'); 
plot(H_bw_ref, M_bw_ref(:,6,6),'-','color',cols(2,:),'displayname', '$\mu_{zz}^{rr}$ PRPB, BW','linewidth',5);
plot(H_sc, M_sc(:,4,4),'o','color',cols(6,:),'displayname', '$\mu_{xx}^{rr}$, SC','markersize',7,'markerfacecolor',cols(6,:))
plot(H_sc_ref, M_sc_ref(:,4,4),'-','color',cols(6,:),'displayname',' $\mu_{xx}^{rr}$ Ref, SC');
plot(H_sc, M_sc(:,6,6),'o','color',cols(7,:),'displayname', '$\mu_{zz}^{rr}$, SC', 'markersize',7,'markerfacecolor',cols(7,:))
plot(H_sc_ref, M_sc_ref(:,6,6),'-','color',cols(7,:),'displayname',' $\mu_{zz}^{rr}$ Ref, SC');

xlim([0,4.8]);
ylim([0,1.01]);
xlabel('$z/R_h$');
ylabel('$(8\pi\eta R_h^3)\mu^{rr}$');
legend show;
legend boxoff;

%% trans-rot
figure(3);
subplot(1,2,1);
semilogy(H_bw, abs(M_bw(:,2,4)),'bs', 'displayname', '$\mu_{yx}^{tr}$ $w=6$'); hold on;
semilogy(H_bw_ref2x, M_bw_ref2x(:,2,4),'o--','color',cols(6,:),'displayname','$\mu_{yx}^{tr}$ $w=12$','MarkerFaceColor',cols(6,:),'MarkerSize', 6)

semilogy(H_bw_ref, abs(M_bw_ref(:,2,4)),'k-','displayname','$\mu_{yx}^{tr}$ PRPB')
semilogy(H_bw_ref_noimg, abs(M_bw_ref_noimg(:,2,4)),'r--','displayname','$\mu_{yx}^{tr}$ RPB')
title('Bottom wall');
xlabel('$z/R_h$');
ylabel('$(6\pi\eta R_h^2)\mu^{tr}$');
ylim([1e-5,1e-1]);
xlim([0,9.6])
set(gca,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1]);
legend show;
legend boxoff;

subplot(1,2,2);
semilogy(H_sc, abs(M_sc(:,2,4)),'bs', 'displayname', '$\mu_{yx}^{tr}$ $w=6$'); hold on;
semilogy(H_sc_ref, abs(M_sc_ref(:,2,4)),'o--', 'color',cols(6,:),'displayname','$\mu_{yx}^{tr}$ $w=12$','MarkerFaceColor',cols(6,:),'MarkerSize', 6)
title('Slit channel');
xlabel('$z/R_h$');
xlim([0,9.6])
ylim([1e-5,1e-1]);
set(gca,'yticklabels',[]);
legend show;
legend boxoff;
