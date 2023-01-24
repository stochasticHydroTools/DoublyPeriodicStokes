clc;clear

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
x_label_y = -0.25;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 1.5;
yl_in = 1.2;
dx_l_in = 0.5;dx_r_in = 0.12;
dy_b_in = 0.45;dy_t_in = 0.25;
xwidth = 3*(xl_in+dx_l_in)+dx_r_in;
ywidth = 2*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1.25 ;
figbox_linewidth = 1;
ms = 3;

COLOR = {[0.6928 0.1651 0.5645] [0.9883 0.6523 0.2114] [0 0.7 0] 'k'};


dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos = zeros(6,4);
pos(1,:) = [0*(dx_l+xl)+dx_l,1*(dy_b+yl+dy_t)+dy_b,xl,yl];
pos(2,:) = [1*(dx_l+xl)+dx_l,1*(dy_b+yl+dy_t)+dy_b,xl,yl];
pos(3,:) = [2*(dx_l+xl)+dx_l,1*(dy_b+yl+dy_t)+dy_b,xl,yl];
pos(4,:) = [0*(dx_l+xl)+dx_l,0*(dy_b+yl+dy_t)+dy_b,xl,yl];
pos(5,:) = [1*(dx_l+xl)+dx_l,0*(dy_b+yl+dy_t)+dy_b,xl,yl];
pos(6,:) = [2*(dx_l+xl)+dx_l,0*(dy_b+yl+dy_t)+dy_b,xl,yl];

R_dp3 = dlmread('data/rollers/Coef_Data/res_scalars_channel_H3_MB_42_L32_q1.txt');
R_dp2 = dlmread('data/rollers/Coef_Data/res_scalars_channel_H3_MB_162_L32_q1.txt');
R_dp1 = dlmread('data/rollers/Coef_Data/res_scalars_channel_H3_MB_642_L32_q1.txt');

R_dp_mb2 = dlmread('data/rollers/Coef_Data/DP_case_res_scalars_channel_H6a_FCM_L32.txt');

R2 = dlmread('data/rollers/Coef_Data/nDP_case_res_scalars_wall_MB_2562_eig_thresh_1xL.txt');
R2mb = dlmread('data/rollers/Coef_Data/nDP_case_res_scalars_wall_MB_1xL.txt');
R2 = flipud(R2);

R2mb_int = R2;
for i = 1:5
    R2mb_int(:,i+1) = interp1(R2mb(:,1),R2mb(:,i+1),R2(:,1),'pchip');
end

h_dp = R_dp3(:,1);
dH = h_dp(2)-h_dp(1);
cut = 2:length(h_dp);
N = length(h_dp);

II = eye(5);
II(3,3) = -1.0;

eta = 1/(6*pi);

R_dp1_sym = eta*[R_dp1(2:N/2+1,2:end); flipud(R_dp1(2:N/2+1,2:end))*II];
R_dp2_sym = eta*[R_dp2(2:N/2+1,2:end); flipud(R_dp2(2:N/2+1,2:end))*II];
R_dp3_sym = eta*[R_dp3(2:N/2+1,2:end); flipud(R_dp3(2:N/2+1,2:end))*II];

[~,indm] = min(abs(R2(:,1)-3));
R_Wall_sym = [R2(1:indm,2:end); flipud(R2(1:indm,2:end))*II];
Rmb_Wall_sym = [R2mb_int(1:indm,2:end); flipud(R2mb_int(1:indm,2:end))*II];

R2x = [R2(1:indm,1); flip(6-R2(1:indm,1))];

%%%%%%%%%%%%%%%%%%%
dRflag = 1;
if dRflag == 0 %0 or 1
   tit_strs = {'$$X^{tt}$$','$$Y^{tt}$$','$$Y^{tr}$$',...
            '$$X^{rr}$$','$$Y^{rr}$$'}; 
end
%%%%%%%%%%%%%%%%%%%


as = [0.0675 0.1310 0.2436];

Rextrap = NaN*R_dp1_sym;
for j = 1:N
    for kk = 1:5
        R_vals = [R_dp1_sym(j,kk) R_dp2_sym(j,kk) R_dp3_sym(j,kk)];
        Rmb_vals = [R_dp_mb2(j,kk+1) R_dp_mb2(j,kk+1) R_dp_mb2(j,kk+1)];
        Rextrap(j,kk) = interp1(as,R_vals-dRflag*Rmb_vals,0,'linear','extrap')+dRflag*R_dp_mb2(j,kk+1);
    end
end

tit_strs = {'$\Delta X^{tt}$','$\Delta Y^{tt}$','$\Delta Y^{tr}$','$\Delta X^{rr}$','$\Delta Y^{rr}$'};
        
smooth_extrap = R_dp3;        

XTICKS = [1,2,3,4,5];
YTICKS = [0     5    10   15    20   ;...
          0     0.2   0.4  0.6   0.8 ;...
         -0.04 -0.02  0    0.02  0.04;...
         -0.04 -0.02  0    0.02  0.04;...
          0     0.2   0.4  0.6   0.8];
YLIMITS = [-2,20;-0.1,0.7;-0.05,0.05;-0.05,0.05;-0.1,0.7];
YTICKLABELS = {'$0$','$5$','$10$','$15$','$20$'     ;...
               '$0$','$0.2$','$0.4$','$0.6$','$0.8$';...
               '$-0.04$','$-0.02$','$0$','$0.02$','$0.04$';...
               '$-0.04$','$-0.02$','$0$','$0.02$','$0.04$';...
               '$0$','$0.2$','$0.4$','$0.6$','$0.8$'};
for i = 1:5
    subplot('position',pos(i,:))
    hold on
    smth = smooth(Rextrap(:,i),'sgolay',2);
    Rextrap(10:end-9,i) = smth(10:end-9);
    smooth_extrap(:,i+1) = Rextrap(:,i);
    hand(4) = plot(R2x,R_Wall_sym(:,i)-dRflag*Rmb_Wall_sym(:,i),'-','color',COLOR{4},'linewidth',linewidth);
    hand(1) = plot(h_dp,R_dp2_sym(:,i)-dRflag*R_dp_mb2(:,i+1),'-','color',COLOR{1},'linewidth',linewidth);
    hand(2) = plot(h_dp,R_dp1_sym(:,i)-dRflag*R_dp_mb2(:,i+1),'--','color',COLOR{2},'linewidth',linewidth);
    hand(3) = plot(h_dp,smooth_extrap(:,i+1)-dRflag*R_dp_mb2(:,i+1),':','color',COLOR{3},'linewidth',linewidth);
    
    
    xlim([1 5])
    ylim(YLIMITS(i,:))

    XLIM = xlim;
    YLIM = ylim;
    xlab = '$h/R_h$';
    x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
%     ylab = '$\mathrm{P}\left(U_x\right)$';
%     x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
%     y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
%     clear text
%     text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)
    
    title(tit_strs{i},'interpreter','latex');
    box on
    set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',XTICKS,'ytick',YTICKS(i,:),'yticklabels',YTICKLABELS(i,:));
    hold off
end

LEGEND = {'SC, $162$ blobs','SC, $642$ blobs','SC, extrapolation','BW, $2562$ blobs'};
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.4*xl,(dy_b+yl+dy_t)+dy_b+0.5*yl,0.05,0.05],'fontsize',small_fontsize);
leg.ItemTokenSize = [14,1];
leg.NumColumns = 1;
legend boxoff

print('rollers_dR_coef','-dpdf','-painters')
