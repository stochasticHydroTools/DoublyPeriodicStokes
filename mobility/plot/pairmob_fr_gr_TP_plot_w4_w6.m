addpath('../');

%this script generate pair mobility function f(r) and g(r) for TP case.
%and validate through comparing with FCM analytical result (ref:maxey 2001)

transparency=0.1; %it sets the traspanrency of the scatter plots for w=6 curves

% load output data from pair_perp_check_TP.py
mobx_u_w6  = load('perpendicular_mobility_unit_w6.txt');
mobx_u_w4  = load('perpendicular_mobility_unit_w4.txt');

mobx_para  = load('parallel_mobility_unit_w4.txt');
mobx_u_w4_gr=mobx_para-mobx_u_w4;

mobx_para  = load('parallel_mobility_unit_w6.txt');
mobx_u_w6_gr=mobx_para-mobx_u_w6;

% set sim params for unit/non_unit grid spacing cases
eta = 1/4/sqrt(pi);
F = 1;

Rh_u_w6 = 1.5539; %1.5539 for w=6; 1.3437 for w=5; 1.2047 for w=4 
Rh_u_w4 = 1.2047;


mu0_u_w6 = 1 / (6 * pi * eta * Rh_u_w6); %mobility in free space
mu0_u_w4 = 1 / (6 * pi * eta * Rh_u_w4); %mobility in free space

Ls_u = linspace(0.,25.,1000); %tried pair distance in [0, 25]

nTrials = 10; %number of trials for each distance Ls_u
Lbox=200;

%compute the anlaytical mob function f(r) with FCM kernel as reference
sigma_w6=Rh_u_w6*sqrt(2.0)/sqrt(pi); %sigma for FCM Gaussians
sigma_w4=Rh_u_w4*sqrt(2.0)/sqrt(pi);

Ls_u(1)=1e-6; %just to temporary remove singularity at zero

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


% absolute difference in pair mobility function f(r) compare with FCM ref
abs_diff_fr_w6 = (mobx_u_w6' +2.8373*Rh_u_w6/Lbox/(6*pi*eta*Rh_u_w6)- fr_FCM);
abs_diff_fr_w4 = (mobx_u_w4' +2.8373*Rh_u_w4/Lbox/(6*pi*eta*Rh_u_w4)- fr_FCM_w4);

abs_diff_gr_w4 = (mobx_u_w4_gr' - gr_FCM_w4);
abs_diff_gr_w6 = (mobx_u_w6_gr' - gr_FCM);

Ls_u(1)=0; %shift back to zero

%plot
tiledlayout(1,2);

ax1 = nexttile;
l1=plot(ax1,Ls_u'/Rh_u_w4, Ls_u'.*(mobx_u_w4+2.8373*Rh_u_w4/Lbox/(6*pi*eta*Rh_u_w4))...
    *8*pi*eta,'c.','MarkerSize',10,'LineWidth',2);hold on;

l2=plot(ax1,Ls_u'/Rh_u_w6, Ls_u'.*(mobx_u_w6+2.8373*Rh_u_w6/Lbox/(6*pi*eta*Rh_u_w6))...
    *8*pi*eta,'r.','MarkerSize',10,'LineWidth',2);hold on;

l3=plot(ax1,Ls_u/Rh_u_w6, Ls_u.*fr_FCM*8*pi*eta,'k-','LineWidth',3);

l4=plot(ax1,Ls_u'/Rh_u_w4, Ls_u'.*mobx_u_w4_gr...
    *8*pi*eta,'g.','MarkerSize',10,'LineWidth',2);hold on;

l5=plot(ax1,Ls_u'/Rh_u_w6, Ls_u'.*mobx_u_w6_gr...
    *8*pi*eta,'m.','MarkerSize',10,'LineWidth',2);hold on;

l6=plot(ax1,Ls_u/Rh_u_w6, Ls_u.*gr_FCM*8*pi*eta,'k--','LineWidth',3); 

hold off
ax1.FontSize=20;
ax1.LineWidth=2;

legend([l1(1) l2(1) l3(1) l4(1) l5(1) l6(1)],{'$f(r)$, $w=4$','$f(r)$, $w=6$','$f_{FCM}(r)$ reference',...
    '$g(r)$, $w=4$','$g(r)$, $w=6$','$g_{FCM}(r)$ reference'}...
    ,'interpreter','latex','FontSize', 26,'Location','southeast')
legend boxoff;

%title(ax1,'Validate pair mobility functions: $f(r)$ and $g(r)$','interpreter','latex','FontSize', 26);
xlabel(ax1,'$r/R_h$','interpreter','latex','FontSize', 26);
ylabel(ax1,'$8\pi\eta r f(r)$ and $8\pi\eta r g(r)$','interpreter','latex','FontSize', 26);

xlim([0 15])
ylim([0 1.16])

%%%%%%%%%%%%%%%plot 2nd scatterred figure%%%%%%%%%%%%%%%%%%%%%%%
ax2 = nexttile;
dot_size=10;

%ZG: these are now useless, just for showing the lengends;
l1=plot(ax2,Ls_u/Rh_u_w4, Ls_u/Rh_u_w4+1e4,'c.','MarkerSize',10,'LineWidth',2);
hold on;

l3=plot(ax2,Ls_u/Rh_u_w4, Ls_u/Rh_u_w4+1e4,'g.','MarkerSize',10,'LineWidth',2);
hold on;

l2=plot(ax2,Ls_u/Rh_u_w6, Ls_u/Rh_u_w6+1e4,'r.','MarkerSize',10,'LineWidth',2);
hold on;

l4=plot(ax2,Ls_u/Rh_u_w6, Ls_u/Rh_u_w6+1e4,'m.','MarkerSize',10,'LineWidth',2);
hold on;

for i=1:10
scatter(Ls_u/Rh_u_w4, abs_diff_fr_w4(i,:),dot_size,...
    'MarkerFaceColor','c','MarkerEdgeColor','c',...
    'MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
hold on;
end

for i=1:10
scatter(Ls_u/Rh_u_w4, abs_diff_gr_w4(i,:),dot_size,...
    'MarkerFaceColor','g','MarkerEdgeColor','g',...
    'MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
hold on;
end

for i=1:10
scatter(Ls_u/Rh_u_w6, abs_diff_fr_w6(i,:),dot_size,...
    'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',transparency,'MarkerEdgeAlpha',transparency);
hold on;
end

for i=1:10
scatter(Ls_u/Rh_u_w6, abs_diff_gr_w6(i,:),dot_size,...
    'MarkerFaceColor','m','MarkerEdgeColor','m',...
    'MarkerFaceAlpha',transparency,'MarkerEdgeAlpha',transparency);
hold on;
end

hold off
ax2.FontSize=20;
ax2.LineWidth=2;
legend([l1(1) l2(1) l3(1) l4(1)],{'$f(r)$, $w=4$','$f(r)$, $w=6$',...
    '$g(r)$, $w=4$','$g(r)$, $w=6$'}...
    ,'interpreter','latex','FontSize', 26,'Location','northeast')
%title(ax2,'Absolute error in $f(r)$ and $g(r)$','interpreter','latex','FontSize', 26);
xlabel(ax2,'$r/R_h$','interpreter','latex','FontSize', 26);
ylabel(ax2,'$f(r)-f_{FCM}(r)$ and $g(r)-g_{FCM}(r)$','interpreter','latex','FontSize', 26);

xlim([0 15])
ylim([-1.5e-3 2.2e-3])

set(gcf,'unit','normalized','position',[0,0,1,0.7]);
legend boxoff;

