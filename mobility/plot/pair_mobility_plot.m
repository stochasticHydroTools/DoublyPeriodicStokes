%% plot pair mobility in either bottom wall or slit channel
%% end of script plots pair mobility matrix asymmetry and posdef-ness

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
addpath('./subtightplot');

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

cols =  distinguishable_colors(20);
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0.1 0.05], [0.1 0.05]);
%% trans-trans
figure(1);
subplot(3,2,1);
pl(1)=plot(H_bw, M_bw(1,:,4,1), 's','color',cols(1,:), 'displayname','FCM $w=6$, $d=3R_h$'); hold on;
pl(2)=plot(H_bw_w4, M_bw_w4(1,:,4,1),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=3R_h$');
pl(3)=plot(H_bw_ref, M_bw_ref(1,:,4,1),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(4)=plot(H_bw, M_bw(2,:,4,1), 's','color',cols(2,:), 'displayname','FCM $w=6$, $d=4R_h$'); hold on;
pl(5)=plot(H_bw_w4, M_bw_w4(2,:,4,1),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=4R_h$');
pl(6)=plot(H_bw_ref, M_bw_ref(2,:,4,1),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(7)=plot(H_bw, M_bw(3,:,4,1), 's','color',cols(6,:), 'displayname','FCM $w=6$, $d=8R_h$'); hold on;
pl(8)=plot(H_bw_w4, M_bw_w4(3,:,4,1),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=8R_h$');
pl(9)=plot(H_bw_ref, M_bw_ref(3,:,4,1),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:9
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h)\mu_{xx}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.44]);
title('Bottom wall');


subplot(3,2,2);
pl(1)=plot(H_sc, M_sc(1,:,4,1), 's','color',cols(1,:), 'displayname','FCM $w=6$, $d=3R_h$'); hold on;
pl(2)=plot(H_sc_w4, M_sc_w4(1,:,4,1),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=3R_h$');
pl(3)=plot(H_sc_ref, M_sc_ref(1,:,4,1),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(4)=plot(H_sc, M_sc(2,:,4,1), 's','color',cols(2,:), 'displayname','FCM $w=6$, $d=4R_h$'); hold on;
pl(5)=plot(H_sc_w4, M_sc_w4(2,:,4,1),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=4R_h$');
pl(6)=plot(H_sc_ref, M_sc_ref(2,:,4,1),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(7)=plot(H_sc, M_sc(3,:,4,1), 's','color',cols(6,:), 'displayname','FCM $w=6$, $d=8R_h$'); hold on;
pl(8)=plot(H_sc_w4, M_sc_w4(3,:,4,1),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=8R_h$');
pl(9)=plot(H_sc_ref, M_sc_ref(3,:,4,1),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:9
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
%ylabel('$\mu_{xx}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.44]);
title('Slit channel');
set(gca,'yticklabel',[])


subplot(3,2,3);
pl(1)=plot(H_bw, M_bw(1,:,6,3), 's','color',cols(1,:), 'displayname','FCM $w=6$, $d=3R_h$'); hold on;
pl(2)=plot(H_bw_w4, M_bw_w4(1,:,6,3),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=3R_h$');
plot(H_bw_ref, M_bw_ref(1,:,6,3),'-','color',cols(1,:),'displayname','Ref, $d=3R_h$');
%
pl(3)=plot(H_bw, M_bw(2,:,6,3), 's','color',cols(2,:), 'displayname','FCM $w=6$, $d=4R_h$'); hold on;
pl(4)=plot(H_bw_w4, M_bw_w4(2,:,6,3),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=4R_h$');
plot(H_bw_ref, M_bw_ref(2,:,6,3),'-','color',cols(2,:),'displayname','Ref, $d=4R_h$');
%
pl(5)=plot(H_bw, M_bw(3,:,6,3), 's','color',cols(6,:), 'displayname','FCM $w=6$, $d=8R_h$'); hold on;
pl(6)=plot(H_bw_w4, M_bw_w4(3,:,6,3),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=8R_h$');
plot(H_bw_ref, M_bw_ref(3,:,6,3),'-','color',cols(6,:),'displayname','Ref, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h)\mu_{zz}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.17]);
set(gca,'ytick',[0,0.05,0.1,0.15]);

subplot(3,2,4);
pl(1)=plot(H_sc, M_sc(1,:,6,3), 's','color',cols(1,:), 'displayname','FCM $w=6$, $d=3R_h$'); hold on;
pl(2)=plot(H_sc_w4, M_sc_w4(1,:,6,3),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=3R_h$');
pl(3)=plot(H_sc_ref, M_sc_ref(1,:,6,3),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(4)=plot(H_sc, M_sc(2,:,6,3), 's','color',cols(2,:), 'displayname','FCM $w=6$, $d=4R_h$'); hold on;
pl(5)=plot(H_sc_w4, M_sc_w4(2,:,6,3),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=4R_h$');
pl(6)=plot(H_sc_ref, M_sc_ref(2,:,6,3),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(7)=plot(H_sc, M_sc(3,:,6,3), 's','color',cols(6,:), 'displayname','FCM $w=6$, $d=8R_h$'); hold on;
pl(8)=plot(H_sc_w4, M_sc_w4(3,:,6,3),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','FCM $w=4$, $d=8R_h$');
pl(9)=plot(H_sc_ref, M_sc_ref(3,:,6,3),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:9
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
%ylabel('$\mu_{zz}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.17]);
set(gca,'yticklabel',[])


subplot(3,2,5);
plot(H_bw, M_bw(1,:,4,3), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(1)=plot(H_bw_w4, M_bw_w4(1,:,4,3),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','$w=4$, $d=3R_h$');
pl(2)=plot(H_bw_ref, M_bw_ref(1,:,4,3),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
plot(H_bw, M_bw(2,:,4,3), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(3)=plot(H_bw_w4, M_bw_w4(2,:,4,3),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','$w=4$, $d=4R_h$');
pl(4)=plot(H_bw_ref, M_bw_ref(2,:,4,3),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
plot(H_bw, M_bw(3,:,4,3), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(5)=plot(H_bw_w4, M_bw_w4(3,:,4,3),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','$w=4$, $d=4R_h$');
pl(6)=plot(H_bw_ref, M_bw_ref(3,:,4,3),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h)\mu_{xz}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.06]);
set(gca,'ytick',[0,0.02,0.04,0.06]);
xlabel('$z/R_h$');

subplot(3,2,6);
%
pl(1)=plot(H_sc, M_sc(1,:,4,3), 's','color',cols(1,:), 'displayname','FCM $w=6$, $d=3R_h$'); hold on;
plot(H_sc_w4, M_sc_w4(1,:,4,3),'s', 'color',cols(1,:),'MarkerFaceColor',cols(1,:), 'MarkerSize',6,'displayname','$w=4$, $d=3R_h$');
pl(2)=plot(H_sc_ref, M_sc_ref(1,:,4,3),'-','color',cols(1,:),'displayname','Reference, $d=3R_h$');
%
pl(3)=plot(H_sc, M_sc(2,:,4,3), 's','color',cols(2,:), 'displayname','FCM $w=6$, $d=4R_h$'); hold on;
plot(H_sc_w4, M_sc_w4(2,:,4,3),'s', 'color',cols(2,:),'MarkerFaceColor',cols(2,:), 'MarkerSize',6,'displayname','$w=4$, $d=4R_h$');
pl(4)=plot(H_sc_ref, M_sc_ref(2,:,4,3),'-','color',cols(2,:),'displayname','Reference, $d=4R_h$');
%
pl(5)=plot(H_sc, M_sc(3,:,4,3), 's','color',cols(6,:), 'displayname','FCM $w=6$, $d=8R_h$'); hold on;
plot(H_sc_w4, M_sc_w4(3,:,4,3),'s', 'color',cols(6,:),'MarkerFaceColor',cols(6,:), 'MarkerSize',6,'displayname','$w=4$, $d=8R_h$');
pl(6)=plot(H_sc_ref, M_sc_ref(3,:,4,3),'-','color',cols(6,:),'displayname','Reference, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
%ylabel('$\mu_{xz}^{tt}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([0,0.06]);
xlabel('$z/R_h$');
set(gca,'yticklabel',[])


for j =1:6
    subplot(3,2,j);
    if j ~= 5 && j ~= 6
        set(gca,'XTick',[]);
    end
end

%% rot-rot
figure(2);

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.01], [0.1 0.05], [0.1 0.05]);
subplot(4,2,1);
%
pl(1)=semilogy(H_bw, M_bw(1,:,11,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,11,8),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,11,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,11,8),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,11,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,11,8),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(8\pi\eta R_h^3)\mu_{yy}^{rr}$','fontsize',33);
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Bottom wall','fontsize',33);
xlim([0,9.6]);
ylim([1e-5,0.03]);

subplot(4,2,3);
%
pl(1)=semilogy(H_bw, M_bw(1,:,10,7), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,10,7),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,10,7), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,10,7),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,10,7), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,10,7),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(8\pi\eta R_h^3)\mu_{xx}^{rr}$','fontsize',33);
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,1e-1]);
set(gca,'ytick',[1e-5, 1e-3])

subplot(4,2,5);
%
pl(1)=semilogy(H_bw, M_bw(1,:,10,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,10,9),'-','color',cols(1,:),'displayname','Ref, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,10,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,10,9),'-','color',cols(2,:),'displayname','Ref, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,10,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,10,9),'-','color',cols(6,:),'displayname','Ref, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(8\pi\eta R_h^3)\mu_{xz}^{rr}$','fontsize',33);
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.015]);


subplot(4,2,7);
%
pl(1)=semilogy(H_bw, M_bw(1,:,12,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,12,9),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,12,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,12,9),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,12,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,12,9),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(8\pi \eta R_h^3)\mu_{zz}^{rr}$','fontsize',33);
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,0.025]);
xlabel('$z/R_h$','fontsize',33);
set(gca,'ytick',[1e-5, 1e-3])


subplot(4,2,2);
%
pl(1)=semilogy(H_sc, M_sc(1,:,11,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref, M_sc_ref(1,:,11,8),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc, M_sc(2,:,11,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref, M_sc_ref(2,:,11,8),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc, M_sc(3,:,11,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref, M_sc_ref(3,:,11,8),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
title('Slit channel','fontsize',33);
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.03]);
set(gca,'yticklabel',[])

subplot(4,2,4);
%
pl(1)=semilogy(H_sc, M_sc(1,:,10,7), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref, M_sc_ref(1,:,10,7),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc, M_sc(2,:,10,7), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref, M_sc_ref(2,:,10,7),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc, M_sc(3,:,10,7), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref, M_sc_ref(3,:,10,7),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,1e-1]);
set(gca,'yticklabel',[])


subplot(4,2,6);
%
semilogy(H_sc(H_sc<=9), M_sc(1,H_sc<=9,10,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
semilogy(H_sc(H_sc<=9), M_sc(2,H_sc<=9,10,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
semilogy(H_sc(H_sc<=9), M_sc(3,H_sc<=9,10,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(1,H_sc_ref<=9,10,9),'-','color',cols(1,:),'displayname','Ref, $d=3R_h$');
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(2,H_sc_ref<=9,10,9),'-','color',cols(2,:),'displayname','Ref, $d=4R_h$');
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(3,H_sc_ref<=9,10,9),'-','color',cols(6,:),'displayname','Ref, $d=8R_h$');
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.015]);
set(gca,'yticklabel',[])

subplot(4,2,8);
%
pl(1)=semilogy(H_sc, M_sc(1,:,12,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref, M_sc_ref(1,:,12,9),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc, M_sc(2,:,12,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref, M_sc_ref(2,:,12,9),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc, M_sc(3,:,12,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref, M_sc_ref(3,:,12,9),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,0.025]);
set(gca,'yticklabel',[])

xlabel('$z/R_h$','fontsize',33);
for j =1:8
    subplot(4,2,j);
    if j ~= 7 && j ~= 8
        set(gca,'XTick',[]);
    end
end


%% trans-rot
figure(3);
subplot(4,2,1);
%
pl(1)=semilogy(H_bw, M_bw(1,:,5,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,5,9),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,5,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,5,9),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,5,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,5,9),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h^2)\mu_{yz}^{tr}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Bottom wall');
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.09]);

subplot(4,2,3);
%
pl(1)=semilogy(H_bw, M_bw(1,:,5,7), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,5,7),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
semilogy(H_bw_ref_noimg, M_bw_ref_noimg(1,:,5,7),'--','color',cols(7,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,5,7), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,5,7),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
semilogy(H_bw_ref_noimg, M_bw_ref_noimg(2,:,5,7),'--','color',cols(8,:),'displayname','RPB, $d=4R_h$');

%
pl(5)=semilogy(H_bw, M_bw(3,:,5,7), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,5,7),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
semilogy(H_bw_ref_noimg, M_bw_ref_noimg(3,:,5,7),'--','color',cols(9,:),'displayname','RPB, $d=8R_h$');

for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h^2)\mu_{yx}^{tr}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,0.01]);

subplot(4,2,5);
%
pl(1)=semilogy(H_bw, M_bw(1,:,4,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,4,8),'-','color',cols(1,:),'displayname','Ref, $d=3R_h$');
pl(3)=semilogy(H_bw_ref_noimg, M_bw_ref_noimg(1,:,4,8),'--','color',cols(7,:),'displayname','RPB, $d=3R_h$');

%
pl(4)=semilogy(H_bw, M_bw(2,:,4,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(5)=semilogy(H_bw_ref, M_bw_ref(2,:,4,8),'-','color',cols(2,:),'displayname','Ref, $d=4R_h$');
pl(6)=semilogy(H_bw_ref_noimg, M_bw_ref_noimg(2,:,4,8),'--','color',cols(8,:),'displayname','RPB, $d=4R_h$');

%
pl(7)=semilogy(H_bw, M_bw(3,:,4,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(8)=semilogy(H_bw_ref, M_bw_ref(3,:,4,8),'-','color',cols(6,:),'displayname','Ref, $d=8R_h$');
pl(9)=semilogy(H_bw_ref_noimg, M_bw_ref_noimg(3,:,4,8),'--','color',cols(9,:),'displayname','RPB, $d=8R_h$');

for j = 1:9
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h^2)\mu_{xy}^{tr}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([6e-6,0.06]);
set(gca,'ytick',[1e-5, 1e-3])

subplot(4,2,7);
%
pl(1)=semilogy(H_bw, M_bw(1,:,6,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_bw_ref, M_bw_ref(1,:,6,8),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_bw, M_bw(2,:,6,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_bw_ref, M_bw_ref(2,:,6,8),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_bw, M_bw(3,:,6,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_bw_ref, M_bw_ref(3,:,6,8),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$(6\pi\eta R_h^2)\mu_{zy}^{tr}$');
hYLabel = get(gca,'YLabel');
%set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.09]);

xlabel('$z/R_h$');


subplot(4,2,2);
%
pl(1)=semilogy(H_sc, M_sc(1,:,5,9), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref, M_sc_ref(1,:,5,9),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc, M_sc(2,:,5,9), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref, M_sc_ref(2,:,5,9),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc, M_sc(3,:,5,9), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref, M_sc_ref(3,:,5,9),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
title('Slit channel');
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.09]);
set(gca,'yticklabel',[])

subplot(4,2,4);
%
pl(1)=semilogy(H_sc(H_sc<=9), M_sc(1,H_sc<=9,5,7), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(1,H_sc_ref<=9,5,7),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc(H_sc<=9), M_sc(2,H_sc<=9,5,7), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(2,H_sc_ref<=9,5,7),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc(H_sc<=9), M_sc(3,H_sc<=9,5,7), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(3,H_sc_ref<=9,5,7),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-6,0.01]);
set(gca,'yticklabel',[])

subplot(4,2,6);
%
semilogy(H_sc(H_sc<=9), M_sc(1,H_sc<=9,4,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
semilogy(H_sc(H_sc<=9), M_sc(2,H_sc<=9,4,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
semilogy(H_sc(H_sc<=9), M_sc(3,H_sc<=9,4,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(1,H_sc_ref<=9,4,8),'-','color',cols(1,:),'displayname','Ref, $d=3R_h$');
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(2,H_sc_ref<=9,4,8),'-','color',cols(2,:),'displayname','Ref, $d=4R_h$');
semilogy(H_sc_ref(H_sc_ref<=9), M_sc_ref(3,H_sc_ref<=9,4,8),'-','color',cols(6,:),'displayname','Ref, $d=8R_h$');
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([6e-6,0.06]);
set(gca,'yticklabel',[])
set(gca,'ytick',[1e-5, 1e-3])
subplot(4,2,8);
%
pl(1)=semilogy(H_sc, M_sc(1,:,6,8), 's','color',cols(1,:), 'displayname','$w=6$, $d=3R_h$'); hold on;
pl(2)=semilogy(H_sc_ref, M_sc_ref(1,:,6,8),'-','color',cols(1,:),'displayname','RPB, $d=3R_h$');
%
pl(3)=semilogy(H_sc, M_sc(2,:,6,8), 's','color',cols(2,:), 'displayname','$w=6$, $d=4R_h$'); hold on;
pl(4)=semilogy(H_sc_ref, M_sc_ref(2,:,6,8),'-','color',cols(2,:),'displayname','RPB, $d=4R_h$');
%
pl(5)=semilogy(H_sc, M_sc(3,:,6,8), 's','color',cols(6,:), 'displayname','$w=6$, $d=8R_h$'); hold on;
pl(6)=semilogy(H_sc_ref, M_sc_ref(3,:,6,8),'-','color',cols(6,:),'displayname','RPB, $d=8R_h$');
for j = 1:6
  set(get(get(pl(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend show;
legend boxoff;
xlim([0,9.6]);
ylim([1e-5,0.09]);
set(gca,'yticklabel',[])
xlabel('$z/R_h$');

for j =1:8
    subplot(4,2,j);
    if j ~= 7 && j ~= 8
        set(gca,'XTick',[]);
    end
end

%% mobility asymmetry and positive definiteness
figure(4);
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.1 0.05], [0.1 0.05]);
set(groot,'defaultLegendFontSize',30);
subplot(2,2,1);
semilogy(fliplr(sort(Masym_bw_w4)),'rs','displayname','BW'); hold on;
semilogy(fliplr(sort(Masym_sc_w4)),'bs','displayname','SC');
title('$w=4$');
ylabel('$\frac{||\mathcal{M}-\mathcal{M}^{T}||_F}{||\mathcal{M}||_F}$');
legend show; legend boxoff;
ylim([1e-13,1e-4]);
set(gca,'XTick',[]);
set(gca,'ytick',[1e-13,1e-10,1e-7,1e-4]);


subplot(2,2,2);
semilogy(fliplr(sort(Masym_bw_w6)),'rs','displayname','BW'); hold on;
semilogy(fliplr(sort(Masym_sc_w6)),'bs','displayname','SC');
title('$w=6$');
set(gca,'XTick',[]);
set(gca,'yticklabel',[]);
ylim([1e-13,1e-4]);

subplot(2,2,3);
semilogy(fliplr(sort(Posdef_bw_w4)), 'rs', 'displayname', 'BW'); hold on;
semilogy(fliplr(sort(Posdef_sc_w4)), 'bs', 'displayname', 'SC');
ylabel('Minimum Eigenvalue of $\mathcal{M}$');
ylim([1e-4,1e0]);

subplot(2,2,4);
semilogy(fliplr(sort(Posdef_bw_w6)), 'rs', 'displayname', 'BW'); hold on;
semilogy(fliplr(sort(Posdef_sc_w6)), 'bs', 'displayname', 'SC');
set(gca,'yticklabel',[]);
%set(gca,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);
ylim([1e-4,1e0]);


