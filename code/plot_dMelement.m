function plot_dMelement(data)

NumberImages=data.In;
selectROImax = size(data.roi{1});
ROIsel=1:selectROImax;

Leftxx=[];Leftxy=[];Leftyy=[];
Rgtxx=[];Rgtxy=[];Rgtyy=[];

for inum=2:NumberImages
    F = data.F{inum};
    
    Ma = data.Mtime2{inum};         % eq.(A1)
    Mb = data.Mtime2{inum-1};

    tmR=data.R{inum};
    tmD=data.D{inum};
    tmA=data.A{inum};
    roi_a      = data.roi{inum};
    roi_b      = data.roi{inum-1};
    cnnt_a   = data.nearC{inum};
    cnnt_b   = data.nearC{inum-1};
        
    MPIV = data.Mpiv{inum-1};

    for i=1:length(ROIsel)
        sel = ROIsel(i);

        link_a = 0;
        indL = roi_a{sel};
        for j=1:length(indL)
            sel2 = indL(j);
            link = logical(sum(cnnt_a==sel2,2));
            link_a = link_a + sum(link);
        end
        link_b = 0;
        indL = roi_b{sel};
        for j=1:length(indL)
            sel2 = indL(j);
            link = logical(sum(cnnt_b==sel2,2));
            link_b = link_b + sum(link);
        end

        dM2 = Ma{sel}/link_a-Mb{sel}/link_b;           % A4 -> M1

        Leftxx = [Leftxx; dM2(1,1)];
        Leftxy = [Leftxy; dM2(1,2)];
        Leftyy = [Leftyy; dM2(2,2)];

        tmpDT = tmR{sel}+tmD{sel}+tmA{sel};

        tm = F{sel}-eye(2)-tmpDT;
        Msave2 = tm*(Mb{sel}/link_b) + (Mb{sel}/link_b)*tm';

        Rgtxx = [Rgtxx;Msave2(1,1)];
        Rgtxy = [Rgtxy;Msave2(1,2)];
        Rgtyy = [Rgtyy;Msave2(2,2)];
    end
end


figh=figure;
figh.Position = [50 100 1000 350];
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile

i4=[-100,100];

h(1)=scatter(Leftxx(:,1),Rgtxx(:,1),20,'b','DisplayName','Q_{xx}');
hold on
h(2)=plot(i4,i4,':k','LineWidth',2.5,'DisplayName','');
xlim([-20 20]); ylim([-20 20]);    xticks([-20 -10 0 10 20]);    yticks([-20 -10 0 10 20]);
axis('square')
legend(h(1),'Location','southeast') ;
str='$$\rm \Delta \hat{M}^{(1)}$$';
xlabel(str,'Interpreter','latex','FontSize',16)
str='$$\rm (\nabla v-D_T)\hat{M}^{(1)}+\hat{M}^{(1)}(\nabla v-D_T)^T$$';
ylabel(str,'Interpreter','latex','FontSize',16)

nexttile

h(1)=scatter(Leftxy(:,1),Rgtxy(:,1),20,'b','DisplayName','Q_{xy}');
hold on
h(2)=plot(i4,i4,':k','LineWidth',2.5,'DisplayName','');
xlim([-20 20]); ylim([-20 20]);    xticks([-20 -10 0 10 20]);    yticks([-20 -10 0 10 20]);
axis('square')
legend(h(1),'Location','southeast') ;
str='$$\rm \Delta \hat{M}^{(1)}$$';
xlabel(str,'Interpreter','latex','FontSize',16)
str='$$\rm (\nabla v-D_T)\hat{M}^{(1)}+\hat{M}^{(1)}(\nabla v-D_T)^T$$';
ylabel(str,'Interpreter','latex','FontSize',16)

nexttile

h(1)=scatter(Leftyy(:,1),Rgtyy(:,1),20,'b','DisplayName','Q_{yy}');
hold on
h(2)=plot(i4,i4,':k','LineWidth',2.5,'DisplayName','');
xlim([-20 20]); ylim([-20 20]);    xticks([-20 -10 0 10 20]);    yticks([-20 -10 0 10 20]);
axis('square')
legend(h(1),'Location','southeast') ;
str='$$\rm \Delta \hat{M}^{(1)}$$';
xlabel(str,'Interpreter','latex','FontSize',16)
str='$$\rm (\nabla v-D_T)\hat{M}^{(1)}+\hat{M}^{(1)}(\nabla v-D_T)^T$$';
ylabel(str,'Interpreter','latex','FontSize',16)


 str=['Figure/exampleFig_kinematics_M.jpg'];
 saveas(figh,str)
 close(figh);

% 
% figh=figure;
% figh.Position = [50 100 1000 350];
% t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% 
% 
% recTrGS=[]; recTrS=[]; recTrR=[]; recTrD=[]; recTrA=[];
% for inum=2:NumberImages
%     R=data.R{inum};    S=data.S{inum};    A=data.A{inum};    D=data.D{inum};    GS=data.Gstar{inum};
% 
%     tmpS=[];    tmpR=[];    tmpD=[];    tmpA=[];    tmpGS=[];
%     for i=1:length(ROIsel)
%         sel= ROIsel(i);
%         % for R
%         trM = trace(R{sel});        tmpR=[tmpR; trM];
%         % for S
%         trM = trace(S{sel});        tmpS=[tmpS; trM];
%         % for GS
%         trM = trace(GS{sel});        tmpGS=[tmpGS; trM];
%         % for D
%         trM = trace(D{sel});        tmpD=[tmpD; trM];
%         % for A
%         trM = trace(A{sel});        tmpA=[tmpA; trM];
%     end
%     recTrR=[recTrR; mean(tmpR)];    recTrS=[recTrS; mean(tmpS)];    recTrGS=[recTrGS; mean(tmpGS)]; recTrD=[recTrD; mean(tmpD)];    recTrA=[recTrA; mean(tmpA)];
% end
% 
% val_MA=12.0;
% ii=1:NumberImages-1;
% MGS = movmean(recTrGS,val_MA);MS = movmean(recTrS,val_MA);MR = movmean(recTrR,val_MA);MD = movmean(recTrD,val_MA);MA = movmean(recTrA,val_MA);
% plot(ii,MGS(:,1),'k','LineWidth',4,'DisplayName','G*')
% hold on
% plot(ii,MS(:,1),'','LineWidth',2,'DisplayName','S')
% plot(ii,MR(:,1),'','LineWidth',2,'DisplayName','R')
% plot(ii,MD(:,1),'','LineWidth',2,'DisplayName','D')
% plot(ii,MA(:,1),'','LineWidth',2,'DisplayName','A')
% xlim([1 NumberImages])
% legend
% strtitle = [num2str(val_MA), ' frames Moving Average'];
% title(strtitle,'FontSize',20)
% xlabel('time','FontSize',20)
% ylabel('Tr(Q)','FontSize',20)
% 
% 
% 
% nexttile
% 
% rhoList=nan(NumberImages,length(ROIsel));
% meanDIV=nan(NumberImages,length(ROIsel));
% trGtmp=nan(NumberImages,length(ROIsel));
% for inum=1:NumberImages-1
%     roi_a      = data.roi{inum};            cnnt_a   = data.nearC{inum};
%     Ma = data.Mtime2{inum}; 
%     MPIV= data.Mpiv{inum};
%     GS=data.Gstar{inum+1};
%     for i=1:length(ROIsel)
%         sel= ROIsel(i);
% 
%         link_a = 0;
%         indL = roi_a{sel};
%         for j=1:length(indL)
%             sel2= indL(j);
%             link = logical(sum(cnnt_a==sel2,2));
%             link_a =link_a + sum(link);
%         end
%         M1=Ma{sel}/link_a;                             
%         rhoList(inum,i) = 1.0/(pi*sqrt(det(M1)));
%         meanDIV(inum,i) = trace(MPIV{sel});
%         trGtmp(inum,i) = trace(GS{sel});
%     end
% end
% 
% XX = ( rhoList(2:NumberImages,:) - rhoList(1:NumberImages-1,:) )./rhoList(1:NumberImages-1,:) + meanDIV(1:NumberImages-1,:);
% YY = ( rhoList(2:NumberImages,:) - rhoList(1:NumberImages-1,:) )./rhoList(1:NumberImages-1,:) + trGtmp(1:NumberImages-1,:);
% 
% meanXX =mean(XX,2,'omitnan');
% meanYY =mean(YY,2,'omitnan');
% 
% trDT=nan(NumberImages-1,length(ROIsel));
% for inum=2:NumberImages
%     tmR=data.R{inum};    tmD=data.D{inum};    tmA=data.A{inum};
%     for i=1:length(ROIsel)
%         sel= ROIsel(i);
%         tmpDT=tmR{sel}+tmD{sel}+tmA{sel};
%         trDT(inum-1,i)=trace(tmpDT);
%     end    
% end
% XX = trDT.*rhoList(1:NumberImages-1,:);
% meanOriTrDT =mean(XX,2,'omitnan');
% meanTrDT =mean(trDT,2,'omitnan');
% 
% 
% ii=1:NumberImages-1;
% 
% MR = movmean(meanXX,val_MA);
% MR2 = movmean(meanYY,val_MA);
% mrTrD = movmean(meanTrDT,val_MA);
% 
% plot(ii,mrTrD,'r','LineWidth',4,'DisplayName','TrDT')
% hold on
% i3=zeros(1,NumberImages-1);
% plot(ii,MR,'b','LineWidth',2,'DisplayName','PIV')
% plot(ii,MR2,'color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'DisplayName','Tr(G)')
% plot(ii,i3,':k','LineWidth',2,'DisplayName','')
% 
% legend
% xlabel('time','FontSize',20)
% ylabel('Tr(D_T)','FontSize',20)
% 
% xlim([1 NumberImages])
% 
% 
% str=['Figure/exampleFig_TrQ_and_TrDT.jpg'];
% saveas(figh,str)
% close(figh);


end