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
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.02], [0.1 0.05], [0.1 0.05]);


pTypes = ["Monopole","Dipole"];

for i = 1:2
    pType = pTypes(i)
    if strcmp(pType,'Monopole')
        load effRad_monopole.mat
    elseif strcmp(pType,'Dipole')
        load effRad_dipole.mat
    end
    Ls = double(Ls); ns = double(ns); ws = double(ws);
    h = Ls(1)./ns(1);
    
    for j = 1:length(ws)
        for k = 1:length(betas)
            for l = 1:nTrials
                if strcmp(pType,'Monopole')
                    p = polyfit(1./Ls',M(:,j,k,l),1);
                    effRh(j,k,l) = 1/h/(6*pi*eta*p(2));
                else
                    p = polyfit(1./Ls'.^3,M(:,j,k,l),1);
                    effRh(j,k,l) = 1/h/(8*pi*eta*p(2))^(1/3);
                end
            end
        end
    end
    Rhs = mean(effRh,3);
    err = 4*std(effRh,0,3)./Rhs*100;
    
    % plot c(w) = R_h/h for each beta
    figure(1);
    subplot(1,2,i);
    cols =  distinguishable_colors(19);
    mspec = ['o','s','d','^','<','v','>','+','h','p','x','*','o','s','.'];
    for k = 1:length(betas)
        spln = spline(ws,Rhs(:,k),linspace(ws(1),ws(end)));
        pl=plot(linspace(ws(1),ws(end)),spln,'-'); hold on;
        pl.Color = cols(k,:);
        set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        for l = 1:nTrials
            pl=plot(ws,effRh(:,k,l),'s','Displayname',strcat('$\beta = $',num2str(betas(k),'%.1f'),'$w$')); grid on;
            pl.Color = cols(k,:);
            if l > 1
                set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end
    xlabel('$w$');
    if i == 1    
        ylabel('$R_h/h=c(w)$');
    else
        legend show; legend boxoff;
        set(gca,'yticklabels',[]);
    end
    title(pType);
    
    % plot c(beta) = R_h/(wh) for each w
    figure(2);
    subplot(1,2,i);
    for j = 1:length(ws)
        pl = plot(betas*ws(j),Rhs(j,:)/(h*ws(j)),'o','Displayname',strcat('$w = $', num2str(ws(j)))); hold on;
        ax = gca; ax.LineWidth = 3;
        pl.Color = cols(j,:);
    end
    %p = polyfit(betas*ws(end),Rhs(end,:)/(h*ws(end)),10);
    betaw = ws'*betas;
    Rhdhw = Rhs./(h*ws');
    p = polyfit(betaw(:),Rhdhw(:),10);
    disp(p');
    bax = linspace(min(betas)*ws(1),max(betas)*ws(end));
    plot(bax, polyval(p,bax), 'k-', 'Displayname', 'polyfit to $c(\beta)$');%strcat('$c(\beta) = $',num2str(p(1)),'$\beta^2 + $', num2str(p(2)), '$\beta + $', num2str(p(3))));

    xlabel('$\beta/w$');
    if i == 1
        ylabel('$R_h/(wh)=c(\beta)$');
        legend show; legend boxoff;
    else
        set(gca,'yticklabels',[]);
    end
    title(pType);
    
    xlabel('$\beta$');
    
    figure(4);    
    subplot(1,2,i);
    for j = 1:length(ws)
        pl = plot(Rhs(j,:)/(h*ws(j)),betas*ws(j),'o','Displayname',strcat('$w = $', num2str(ws(j)))); hold on;
        pl.Color = cols(j,:);
    end
    %invp = polyfit(polyval(p,union(betas*ws(1),betas*ws(end))),union(betas*ws(1),betas*ws(end)),5);
    invp = polyfit(Rhdhw(:),betaw(:),10);

    disp(invp');
    bax_inv = linspace(Rhs(end,end)/(ws(end)*h), Rhs(1,1)/(ws(1)*h));
    plot(bax_inv, polyval(invp,bax_inv), 'k-', 'Displayname', 'polyfit to $c^{-1}(\beta)$')        
    
    xlabel('$R_h/(wh)=c(\beta)$');
    ylabel('$\beta$');
    title(pType);
    legend show; legend boxoff;
    
    if strcmp(pType,'Monopole') 
        dlmwrite('monopole_cbeta.txt',p','precision',16);
        dlmwrite('monopole_cbeta_inv.txt',invp','precision',16);
    elseif strcmp(pType,'Dipole')
        dlmwrite('dipole_cbeta.txt',p','precision',16);
        dlmwrite('dipole_cbeta_inv.txt',invp','precision',16);
    end  
    
    % plot percent error in c(w)
    figure(3); 
    bax = linspace(min(betas),max(betas));
    subplot(1,2,i);
    for j = 1:length(ws)
        splnE(j) = pchip(betas*ws(j),err(j,:));
        semilogy(betas,err(j,:),'s','color',cols(j,:),'Displayname',strcat('$w = $',num2str(ws(j)))); hold on;
        ax = gca; ax.LineWidth = 3; 
        bax_err = linspace(min(betas*ws(j)),max(betas)*ws(j),1000);
        splneval = ppval(splnE(j),bax_err);
        if strcmp(pType,'Monopole') 
            if j == 1
                dlmwrite('monopole_w4_err_spline.txt',[bax_err',splneval'],'precision',16);
            elseif j == 2
                dlmwrite('monopole_w5_err_spline.txt',[bax_err',splneval'],'precision',16);
            elseif j == 3
                dlmwrite('monopole_w6_err_spline.txt',[bax_err',splneval'],'precision',16);
            end
        else
            if j == 1
                dlmwrite('dipole_w4_err_spline.txt',[bax_err',splneval'],'precision',16);
            elseif j == 2
                dlmwrite('dipole_w5_err_spline.txt',[bax_err',splneval'],'precision',16);
            elseif j == 3
                dlmwrite('dipole_w6_err_spline.txt',[bax_err',splneval'],'precision',16);
            end
        end
        pl=semilogy(bax,ppval(pchip(betas,err(j,:)),bax),'-','color',cols(j,:));
        set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('$\beta/w$');
    if i == 1
        ylabel('\%-error','rotation',90);
        legend show; legend boxoff;
    else
        set(gca,'yticklabels',[]);
    end
    title(pType);
end
