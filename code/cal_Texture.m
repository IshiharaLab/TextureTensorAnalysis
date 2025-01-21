function datas=cal_Texture(handles)

% calculate the M

NumberImages = handles.In;
roi      = handles.roi{1};
d_roi        = size(roi,2);
 

 disp("2.3.1 calculation of Mpiv")

% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  Matrix of PIV (D+\Omega)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recMPIV=cell(NumberImages-1,1);
recDPIV=cell(NumberImages-1,1);

roiT = handles.roi;

se = strel('diamond',1);
x=handles.pivx;
y=handles.pivy;
u=handles.pivu;
v=handles.pivv;


[lx,ly]=size(x{1});
dx = x{1}(1,2)-x{1}(1,1);
for inum = 1:NumberImages-1
  cen_x = x{inum}(2:lx-1,2:ly-1);
  cen_y = y{inum}(2:lx-1,2:ly-1);
  ux = ( u{inum}(2:lx -1, 3:ly) - u{inum}(2:lx -1, 1:ly-2) )/(2.0*dx);
  uy = ( u{inum}(3:lx, 2:ly-1) - u{inum}(1:lx-2, 2:ly-1) )/(2.0*dx);
  vx = ( v{inum}(2:lx -1, 3:ly) - v{inum}(2:lx -1, 1:ly-2) )/(2.0*dx);
  vy = ( v{inum}(3:lx, 2:ly-1) - v{inum}(1:lx-2, 2:ly-1) )/(2.0*dx);

    disp("frame step in PIV :")
    disp(inum)

  MPIV=cell(d_roi,1);
  DPIV=cell(d_roi,1);
  xq= cen_x(:);
  yq= cen_y(:);
  areaA=handles.AreaLabel{inum};
  idROI=roiT{inum};
  for jroi=1:d_roi
      AA =ismember(areaA,idROI{jroi});
      A3=imclose(AA,se);
      blocations = bwboundaries(A3,'noholes');
      postmp = blocations{1};
      ROIedge = fliplr(postmp);
    [in,~] = inpolygon(xq,yq,ROIedge(:,1),ROIedge(:,2));
    MPIV{jroi}=[mean(ux(in)), mean(uy(in)) ; mean(vx(in)), mean(vy(in)) ] ;
    DPIV{jroi}=[mean(ux(in)), 0.5*( mean(uy(in))+mean(vx(in)) ) ; 0.5*( mean(uy(in))+mean(vx(in)) ) , mean(vy(in)) ] ;
  end
  recMPIV{inum}=MPIV;
  recDPIV{inum}=DPIV;
end

handles.Mpiv = recMPIV;
handles.Dpiv = recDPIV;
 

 disp("2.3.2 calculation of M")



recMtime = cell(NumberImages,1);
recMtime2 = cell(NumberImages,1);
recG_mc = cell(NumberImages,1);
recG_Mc = cell(NumberImages,1);

recRind = cell(NumberImages,1);
recrind = cell(NumberImages,1);

recS = cell(NumberImages,1);
recR = cell(NumberImages,1);
recD = cell(NumberImages,1);
recA = cell(NumberImages,1);
recJ = cell(NumberImages,1);
nnS = cell(d_roi,1);
nnR = cell(d_roi,1);
nnD = cell(d_roi,1);
nnA = cell(d_roi,1);
nnJ = cell(d_roi,1);

recmT = cell(NumberImages,1);
recMT = cell(NumberImages,1);
rechLinc = cell(NumberImages,1);        % [#nc #nT #NT]

for jroi=1:d_roi
    nnS{jroi} = [];
    nnR{jroi} = [];
    nnD{jroi} = [];
    nnA{jroi} = [];
    nnJ{jroi} = [];
end



souka=[];
soujou=[];

recF0       = cell(NumberImages,1);
recF1       = cell(NumberImages,1);
recF         = cell(NumberImages,1);
recGstar   = cell(NumberImages,1);
recPsi   = cell(NumberImages,1);


recG_mc{1}   = nan;
recG_Mc{1}   = nan;
recS{1}  = nan;
recR{1}  = nan;
recD{1}  = nan;
recA{1}  = nan;
recJ{1}  = nan;
recF0{1} = nan;
recF1{1} = nan;
recF{1}  = nan;
recGstar{1}  = nan;
recmT{1}  = nan;
recMT{1}  = nan;
rechLinc{1}  = nan;        % [ #nc #nT #NT]

% 1 frame me
inum    = 1;
xq       = handles.pos{inum}(:,2);
yq       = handles.pos{inum}(:,1);
cnnt   = handles.nearC{inum};
Mrec   = cell(d_roi,1);

% vtx4
vt4     = handles.vtx4{inum};
% extended 4vtx
if isfield(handles,'new4vtx')
    a = [vt4;handles.new4vtx{inum}];
    vt4 = unique(a,'rows');
end

for jroi=1:d_roi
    indL = roi{jroi};
    M=zeros(2,2);
    num=0;
    for j = 1:length(indL)
        sel = indL(j);
        [M2,~] = Mtensor(M,num,sel,cnnt,xq,yq,vt4);
        M=M2;
    end
    Mrec{jroi} = M;
end
recMtime{inum} = Mrec;

Mrec2 = cell(d_roi,1);
for jroi=1:d_roi
    indL = roi{jroi};
    M=zeros(2,2);
    num=0;
    for j = 1:length(indL)
        sel = indL(j);
        [M2,~] = Mtensor_old(M,num,sel,cnnt,xq,yq,vt4);
        M=M2;
    end
    Mrec2{jroi} = M;
end
recMtime2{inum} = Mrec2;


for inum = 2:NumberImages
    disp("frame step in strain :")
    disp(inum)
    
    xa      = handles.pos{inum}(:,2);
    ya      = handles.pos{inum}(:,1);
    xb      = handles.pos{inum-1}(:,2);
    yb      = handles.pos{inum-1}(:,1);
    ind1    = handles.ind1{inum-1};
    ind2    = handles.ind2{inum-1};
    cnnt_b  = handles.nearC{inum-1};
    cnnt    = handles.nearC{inum};
    
    vt4b    = handles.vtx4{inum-1};
    vt4a    = handles.vtx4{inum};
    oCell  = handles.outcell{inum};

    if isfield(handles,'new4vtx')
        a = [vt4a;handles.new4vtx{inum}];
        vt4a = unique(a,'rows');
        a = [vt4b;handles.new4vtx{inum-1}];
        vt4b = unique(a,'rows');
    end

    
    
    roi_a      = handles.roi{inum};
    roi_b      = handles.roi{inum-1};
    % D,A,J_edge index  % check
    if isempty(handles.Divtime{inum-1})
        cell_bef_D=nan;
    else
        cell_bef_D = handles.Divtime{inum-1}(:,1);
    end
    
    ind_Rlink=[];
    ind_rlink=[];
    
    NumDcell= cell(d_roi,1);
    
    mTrec= cell(d_roi,1);
    MTrec= cell(d_roi,1);
    numLinkrec= cell(d_roi,1);
    
    Mrec= cell(d_roi,1);
    Mrec2 = cell(d_roi,1);
    G_mc= cell(d_roi,1);
    G_Mc= cell(d_roi,1);
    GS  = cell(d_roi,1);
    tPsi = cell(d_roi,1);
    S   = cell(d_roi,1);
    R   = cell(d_roi,1);
    D   = cell(d_roi,1);
    Ap  = cell(d_roi,1);
    J   = cell(d_roi,1);
    F0  = cell(d_roi,1);
    F1  = cell(d_roi,1);
    F   = cell(d_roi,1);
    for jroi=1:d_roi
        % n series are the number of half-link to correlate (1) before frame,(2) after frame   nJ(1)=nJ(1)+length(p)
        nG = zeros(1,2);
        nR = zeros(1,2);
        nD = zeros(1,2);
        nA = zeros(1,2);
        nJ = zeros(1,2);
        
        %list in jroi-th ROI
        indL = roi_a{jroi};
        indLb = roi_b{jroi};
                
        % h-link sousuu in ROI
        tot_Link = 0;
        for j=1:length(indL)
            sel= indL(j);
            link = logical(sum(cnnt==sel,2));
            tmpM = cnnt(link,:);
            each_ind = tmpM(tmpM~=sel);
            tot_Link = tot_Link + length(each_ind);
        end
        bef_Link = 0;
        for j=1:length(indLb)
            sel= indLb(j);
            link = logical(sum(cnnt_b==sel,2));
            tmpM = cnnt_b(link,:);
            each_ind = tmpM(tmpM~=sel);
            bef_Link = bef_Link + length(each_ind);
        end
        

        % M calculation (all cell in ROI)
        num=0;
        M   = zeros(2,2);
        for j=1:length(indL)
            sel= indL(j);
            [M2,~] = Mtensor(M,num,sel,cnnt,xa,ya,vt4a);
            M=M2;
        end
        Mrec{jroi} = M;

        M   = zeros(2,2);
        for j=1:length(indL)
            sel= indL(j);
            [M2,~] = Mtensor_old(M,num,sel,cnnt,xa,ya,vt4a);
            M=M2;
        end
        Mrec2{jroi} = M;
        
        % apotothis   (lost cell between before and now frame)
        MMA=zeros(2,2);
        numtmp=0;
        tmp = ind1(ismember(ind1,indLb));
        lost_b = indLb(~ismember(indLb,tmp));
        if isempty(lost_b)==0
            for j=1:length(lost_b)
                sel_b= lost_b(j);
                [M2,numtmp2] = Mtensor_old(MMA,numtmp,sel_b,cnnt_b,xb,yb,vt4b);
                MMA=M2;
                numtmp=numtmp2;
            end
            nA(1) = numtmp;
        end
                
        % division
        MmD=zeros(2,2);
        MMD=zeros(2,2);
        D_b = indLb(ismember(indLb,cell_bef_D));
        NumDcell{jroi}=length(D_b);
        if isempty(D_b)==0
            tmp = handles.Divtime{inum-1}(:,:);
            tmp2 = tmp(ismember(cell_bef_D,D_b),2:3);
            D_a = tmp2(:);
            
            numtmp=0;
            for j=1:length(D_b)
                sel_b= D_b(j);
                [M2,numtmp2] = Mtensor_old(MMD,numtmp,sel_b,cnnt_b,xb,yb,vt4b);
                MMD=M2;
                numtmp=numtmp2;
            end
            nD(1) = numtmp;
                
            numtmp=0;
            for j=1:length(D_a)
                sel= D_a(j);
                [M2,numtmp2] = Mtensor_old(MmD,numtmp,sel,cnnt,xa,ya,vt4a);
                MmD=M2;
                numtmp=numtmp2;
            end
            nD(2) = numtmp;
        end
        
        % for J  % ONLY out of region in Lagragian
        MmJ=zeros(2,2);
        MMJ=zeros(2,2);
        tmp = ind2(ismember(ind1,indLb));
        a = tmp(ismember(tmp,oCell));       
        lost_b = ind1(ismember(ind2,a));
        if isempty(lost_b)==0
            numtmp=0;
            for j=1:length(lost_b)
                sel_b= lost_b(j);
                [M2,numtmp2] = Mtensor_old(MMJ,numtmp,sel_b,cnnt_b,xb,yb,vt4b);
                MMJ=M2;
                numtmp=numtmp2;
            end
            nJ(1) = numtmp;
        end
      
        % remove the division cell in indL
        aa = handles.Divtime{inum-1}(:,:);
        if   ~isempty(aa)
            tmp = indL(~ismember(indL,handles.Divtime{inum-1}(:,2:3)));         % eject the division
            indL = tmp;
        end
        
        
        D_list=handles.Divtime{inum-1};
        Mmc = zeros(2,2);
        MMc = zeros(2,2);
        MmR = zeros(2,2);
        MMR = zeros(2,2);
        
        Fll = zeros(2,2);
        FLl = zeros(2,2);
        FlL = zeros(2,2);
        FLL = zeros(2,2);
        ncforrec =0;
        for j=1:length(indL)
            sel= indL(j);
            link=logical(sum(cnnt==sel,2));
            tmpM=cnnt(link,:);
            now_cnt=tmpM(tmpM~=sel);
            bef_nei=ind1(ismember(ind2,now_cnt));
            sel_b=ind1(ind2==sel);                                          % id of before frame of selected cell
            if (  isempty(sel_b)==0 && ismember(sel_b,cell_bef_D)==0  )     % (there are tracking id) \hat (nothing in id of division cells)
                link_b=logical(sum(cnnt_b==sel_b,2));
                tmp=cnnt_b(link_b,:);
                bef_cnt=tmp(ismember(cnnt_b(link_b,:),sel_b)~=1);
                
                lost_cnt = bef_cnt(~ismember(bef_cnt,bef_nei));     % lost link
                keep_cnt = bef_cnt(ismember(bef_cnt,bef_nei));      % Remaining links
                keep_now = [];
                for j4=1:length(keep_cnt)
                    tmp = keep_cnt(j4);
                    keep_now =  [keep_now; ind2(ismember(ind1,tmp)) ];
                end
                new_cnt = now_cnt(~ismember(now_cnt,keep_now));           % new link
                
                         
                if ~isempty(lost_cnt(ismember(lost_cnt,cell_bef_D)))
                    tmp_lost=lost_cnt(~ismember(lost_cnt,cell_bef_D));
                    tmp=lost_cnt(ismember(lost_cnt,cell_bef_D));
                    for jj_D =1:length(tmp)
                        tmp_sel=tmp(jj_D);
                        
                        aa=ismember(D_list(:,1),tmp_sel);
                        kouho=D_list(aa,2:3);
                        kouho_now=now_cnt(ismember(now_cnt,kouho));
                        if length(kouho_now)==0
                            tmp_lost=[tmp_lost;tmp_sel];
                        else  
                            LLx=xb(tmp_sel)-xb(sel_b);         % lost link
                            LLy=yb(tmp_sel)-yb(sel_b);
                            llx=xa(kouho_now)-xa(sel);                % new link
                            lly=ya(kouho_now)-ya(sel);
                            MMD = MMD + 0.5*[sum(LLx.*LLx), sum(LLx.*LLy); sum(LLx.*LLy), sum(LLy.*LLy)];
                            MmD = MmD + 0.5*[sum(llx.*llx), sum(llx.*lly); sum(llx.*lly), sum(lly.*lly)];
                            nD(1) = nD(1) + 1;
                            nD(2) = nD(2) + length(kouho_now);
                            new_cnt = new_cnt(~ismember(new_cnt,kouho_now));
                        end
                    end
                    lost_cnt=tmp_lost;
                end
                       
                
                
                llx=xa(new_cnt)-xa(sel);                % new link
                lly=ya(new_cnt)-ya(sel);
                LLx=xb(lost_cnt)-xb(sel_b);         % lost link
                LLy=yb(lost_cnt)-yb(sel_b);

                if  ~isempty(lost_cnt)
                    for ioo=1:length(lost_cnt)
                        ind_Rlink=[ind_Rlink;lost_cnt(ioo),sel_b];
                    end
                end
                
                if  ~isempty(new_cnt)
                    for ioo=1:length(new_cnt)
                        ind_rlink=[ind_rlink;new_cnt(ioo),sel];
                    end
                end
                
                lx=xa(keep_now)-xa(sel);            % Links left
                ly=ya(keep_now)-ya(sel);
                Lx=xb(keep_cnt)-xb(sel_b);      % Remaining links
                Ly=yb(keep_cnt)-yb(sel_b);

                ncforrec = ncforrec + length(keep_cnt);
                
                % in-in G and F
                if ismember(sel,indL)==1
                    tmp=length(Lx);
                    MMc = MMc + 0.5*[sum(Lx.*Lx), sum(Lx.*Ly); sum(Lx.*Ly), sum(Ly.*Ly)];           % Mc in G
                    nG(1) = nG(1) + tmp;
                    tmp=length(lx);
                    Mmc = Mmc + 0.5*[sum(lx.*lx), sum(lx.*ly); sum(lx.*ly), sum(ly.*ly)];           % mc in G
                    nG(2) = nG(2) + tmp;
                    
                    tmpRN =length(lost_cnt);
                    tmpRn =length(new_cnt);
                    %   re-arrengement index
                    MMR2 = 0.5*[sum(LLx.*LLx), sum(LLx.*LLy); sum(LLx.*LLy), sum(LLy.*LLy)];
                    MmR2 = 0.5*[sum(llx.*llx), sum(llx.*lly); sum(llx.*lly), sum(lly.*lly)];
                    

                    % for 4vtx before
                    chk4 = sum(ismember(vt4b,sel_b),1:2);
                    if chk4 ~= 0
                        a=logical(sum(ismember(vt4b,sel_b),2));
                        tmp4vtx = vt4b(a,:);
                        each4vt=[bef_cnt;sel_b];
                        tmp2 = tmp4vtx(~ismember(tmp4vtx,each4vt));

                        for i_vt4 =1:length(tmp2)
                            tmp3=tmp2(i_vt4);
                            lxtm=xb(tmp3)-xb(sel_b);
                            lytm=yb(tmp3)-yb(sel_b);
                            BB = ind2(ismember(ind1,tmp3));
                            if max(sum(ismember(cnnt,[BB,sel]),2))==2
                                MMc = MMc + 0.25*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                nG(1) = nG(1) + 0.5;

                                lxtm2=xa(BB)-xa(sel);
                                lytm2=ya(BB)-ya(sel);
                                Mmc = Mmc + 0.25*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                nG(2) = nG(2) + 0.5;
                                % F0, F1, F
                                Fll=Fll + 0.5*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                FLl=FLl + 0.5*[lxtm2*lxtm, lxtm*lytm2; lxtm2*lytm, lytm2*lytm];
                                FlL=FlL + 0.5*[lxtm2*lxtm, lxtm2*lytm; lxtm*lytm2, lytm2*lytm];
                                FLL=FLL + 0.5*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];

                                MmR2 = MmR2 - 0.25*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                tmpRn = tmpRn - 0.5;
                            else
                                if max(sum(ismember(vt4a,[BB,sel]),2))==2       % from 4vtx to 4vtx
                                    MMc = MMc + 0.25*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                    nG(1) = nG(1) + 0.5;

                                    lxtm2=xa(BB)-xa(sel);
                                    lytm2=ya(BB)-ya(sel);
                                    Mmc = Mmc + 0.25*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                    nG(2) = nG(2) + 0.5;
                                    % F0, F1, F
                                    Fll=Fll + 0.5*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                    FLl=FLl + 0.5*[lxtm2*lxtm, lxtm*lytm2; lxtm2*lytm, lytm2*lytm];
                                    FlL=FlL + 0.5*[lxtm2*lxtm, lxtm2*lytm; lxtm*lytm2, lytm2*lytm];
                                    FLL=FLL + 0.5*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                else
                                    MMR2 = MMR2+0.25*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                    tmpRN = tmpRN + 0.5;
                                end
                            end
                        end
                    end

                    % for 4vtx after
                    chk4 = sum(ismember(vt4a,sel),1:2);
                    if chk4 ~= 0
                        a=logical(sum(ismember(vt4a,sel),2));
                        tmp4vtx = vt4a(a,:);
                        each4vt=[now_cnt;sel];
                        tmp2 = tmp4vtx(~ismember(tmp4vtx,each4vt));

                        for i_vt4 =1:length(tmp2)
                            tmp3=tmp2(i_vt4);
                            lxtm2=xa(tmp3)-xa(sel);
                            lytm2=ya(tmp3)-ya(sel);
                            AA = ind1(ismember(ind2,tmp3));
                            if max(sum(ismember(cnnt,[AA,sel_b]),2))==2
                                Mmc = Mmc + 0.25*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                nG(2) = nG(2) + 0.5;

                                lxtm=xb(AA)-xb(sel_b);
                                lytm=yb(AA)-yb(sel_b);
                                MMc = MMc + 0.25*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                nG(1) = nG(1) + 0.5;
                                MMR2 = MMR2 - 0.25*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                                tmpRN = tmpRN - 0.5;
                                % F0, F1, F
                                Fll=Fll + 0.5*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                FLl=FLl + 0.5*[lxtm2*lxtm, lxtm*lytm2; lxtm2*lytm, lytm2*lytm];
                                FlL=FlL + 0.5*[lxtm2*lxtm, lxtm2*lytm; lxtm*lytm2, lytm2*lytm];
                                FLL=FLL + 0.5*[lxtm*lxtm, lxtm*lytm; lxtm*lytm, lytm*lytm];
                            else
                                MmR2 = MmR2+0.25*[lxtm2*lxtm2, lxtm2*lytm2; lxtm2*lytm2, lytm2*lytm2];
                                tmpRn = tmpRn + 0.5;
                            end
                        end
                    end
                    
                    
                    if tmpRN ~= 0
                        nR(1) = nR(1) + tmpRN;
                        MMR = MMR + MMR2;        
                    end
                    if tmpRn ~= 0
                        MmR = MmR + MmR2;
                        nR(2) = nR(2) + tmpRn;
                    end
                    
                    % F0, F1, F
                    Fll=Fll + [sum(lx.*lx), sum(lx.*ly); sum(ly.*lx), sum(ly.*ly)];
                    FLl=FLl + [sum(Lx.*lx), sum(Lx.*ly); sum(Ly.*lx), sum(Ly.*ly)];
                    FlL=FlL + [sum(lx.*Lx), sum(lx.*Ly); sum(ly.*Lx), sum(ly.*Ly)];
                    FLL=FLL + [sum(Lx.*Lx), sum(Lx.*Ly); sum(Ly.*Lx), sum(Ly.*Ly)];
                end
            end
        end
        
        numLinkrec{jroi} = [ ncforrec, tot_Link-ncforrec, bef_Link-ncforrec];  % [#nc #nT #NT]
        
        G_mc{jroi}=Mmc;
        G_Mc{jroi}=MMc;
        if det(FLL)==0 || det(FLl)==0
            F0{jroi}=eye(2);
            F1{jroi}=eye(2);
        else
            F0{jroi}=FlL*inv(FLL);
            F1{jroi}=Fll*inv(FLl);
        end
        
        % strenge of tensor A   sqrt( 0.5*sum(A.*A,'all') )
         % F{jroi} = sqrt( sqrt(0.5*(sum(F1{jroi}.*F1{jroi},'all')))/sqrt(0.5*(sum(F0{jroi}.*F0{jroi},'all'))) )*F0{jroi};     %  Geometric mean
        F{jroi}=0.5*(F0{jroi}+F1{jroi});        % Arithmetic mean
        Ftm = F{jroi};
        Gstar2= 0.5*(Ftm'*Ftm-eye(2));          
        %    Gstar = Mmc-MMc-psi;                             
        
        %apporox check
        tmpG=Mmc-MMc;
        tmpbunsi=Mmc-Ftm*MMc*Ftm';
        tmpA1=sqrt( 0.5*(sum(tmpG.*tmpG,'all')) );
        tmpB1=sqrt( 0.5*(sum(tmpbunsi.*tmpbunsi,'all')) );
        souka = [souka;tmpB1/tmpA1];
        
        Ftm2 = sqrt( sqrt(0.5*(sum(F1{jroi}.*F1{jroi},'all')))/sqrt(0.5*(sum(F0{jroi}.*F0{jroi},'all'))) )*F0{jroi};
        %         Ftm2 = sqrt( sqrt(0.5*(sum(F0{jroi}.*F0{jroi},'all')))/sqrt(0.5*(sum(F1{jroi}.*F1{jroi},'all'))) )*F1{jroi};
        tmpbunsi2=Mmc-Ftm2*MMc*Ftm2';
        tmpB2=sqrt( 0.5*(sum(tmpbunsi2.*tmpbunsi2,'all')) );
        soujou = [soujou;tmpB2/tmpA1];
        
        
        GS{jroi} = Gstar2;
        %tilda check
        if sum(MMc(:))==0
            S{jroi}= zeros(2,2);
            R{jroi}= zeros(2,2);
            D{jroi}= zeros(2,2);
            Ap{jroi}= zeros(2,2);
            J{jroi}= zeros(2,2);
            
            mTrec{jroi}=zeros(2,2);
            MTrec{jroi}=zeros(2,2);
        else
            psi = tilda_psi(Ftm,MMc);            
            tPsi{jroi}=psi;

            tmp = Mmc + MmR + MmD + MmJ;
            til_m = tilda_Q(tmp,MMc);
            tmp = MMc + MMR + MMD + MMA + MMJ;
            til_M = tilda_Q(tmp,MMc);
            
            S{jroi} = ( (nG(1)+nR(1)+nD(1)+nA(1)+nJ(1))/(nG(2)+nR(2)+nD(2)+nA(2)+nJ(2)) )*til_m - til_M - psi;                     
            tmp1 = tilda_Q(MmR,MMc);
            tmp2 = tilda_Q(MMR,MMc);
            R{jroi} = ( (nR(2)-nR(1))/(nG(2)+nR(2)+nD(2)+nA(2)+nJ(2)) )*til_m - (tmp1-tmp2);
            
            tmp1 = tilda_Q(MmD,MMc);
            tmp2 = tilda_Q(MMD,MMc);
            D{jroi} = ( (nD(2)-nD(1))/(nG(2)+nR(2)+nD(2)+nA(2)+nJ(2)) )*til_m - (tmp1-tmp2);
            
            tmp2 = tilda_Q(MMA,MMc);
            Ap{jroi} = ( (nA(2)-nA(1))/(nG(2)+nR(2)+nD(2)+nA(2)+nJ(2)) )*til_m + tmp2;
            
            tmp1 = tilda_Q(MmJ,MMc);
            tmp2 = tilda_Q(MMJ,MMc);
            J{jroi} = ( (nJ(2)-nJ(1))/(nG(2)+nR(2)+nD(2)+nA(2)+nJ(2)) )*til_m - (tmp1-tmp2);
            
            mTrec{jroi} = MmR + MmD + MmJ;
            MTrec{jroi} = MMR + MMD + MMA + MMJ;    
        end
        
        nnS{jroi} = [nnS{jroi};nG(1) nG(2)];
        nnR{jroi} = [nnR{jroi};nR(1) nR(2)];
        nnD{jroi} = [nnD{jroi};nD(1) nD(2)];
        nnA{jroi} = [nnA{jroi};nA(1) nA(2)];
        nnJ{jroi} = [nnJ{jroi};nJ(1) nJ(2)];
    end
    recMtime{inum}=Mrec;
    recMtime2{inum}=Mrec2;
    recF0{inum}= F0;
    recF1{inum}= F1;
    recF{inum}= F;
    
    recG_mc{inum}= G_mc;
    recG_Mc{inum}= G_Mc;
    recGstar{inum}= GS;
    recPsi{inum}=tPsi;

    recS{inum}= S;
    recR{inum}= R;
    recD{inum}= D;
    recA{inum}= Ap;
    recJ{inum}= J;
    
    recmT{inum} = mTrec;
    recMT{inum} = MTrec;
    
    recRind{inum} = ind_Rlink;
    recrind{inum} = ind_rlink;
    rechLinc{inum} = numLinkrec;

    recNumDcell{inum} = NumDcell;
    
    cla
    mask_sepa=handles.IROI(:,:,inum);
    RGB2 = label2rgb(mask_sepa,'pink','c','shuffle');
    clear mask_sepa;
    imshow(RGB2);
    hold on;
    
    xcent = handles.roipos{inum};
    tmpt=0:10:360;
    for jroi=1:d_roi
        [Vec,Dval] = eig(Mrec{jroi});
        Dval=0.05*Dval;
        x1=[sqrt(Dval(1,1))*cosd(tmpt);sqrt(Dval(2,2))*sind(tmpt)];
        ctmp=[xcent{jroi}(1)*ones(1,length(tmpt));xcent{jroi}(2)*ones(1,length(tmpt))];
        X=Vec*x1+ctmp;
        plot(X(1,:),X(2,:),"r","LineWidth",2.5);
        
        lng=[xcent{jroi}(1) xcent{jroi}(2); xcent{jroi}(1)  xcent{jroi}(2)] +[sqrt(Dval(1,1))*Vec(1,:); -sqrt(Dval(1,1))*Vec(1,:)];
        lng2=[xcent{jroi}(1) xcent{jroi}(2); xcent{jroi}(1)  xcent{jroi}(2)] +[sqrt(Dval(2,2))*Vec(2,:); -sqrt(Dval(2,2))*Vec(2,:)];
        if  Dval(1,1) > Dval(2,2)
            plot(lng(:,1),lng(:,2),'r',"LineWidth",1.5);
            plot(lng2(:,1),lng2(:,2),'b',"LineWidth",1.5);
        else
            plot(lng(:,1),lng(:,2),'b',"LineWidth",1.5);
            plot(lng2(:,1),lng2(:,2),'r',"LineWidth",1.5);
        end
    end
    
    hold off;
end

handles.indR=recRind;
handles.indr=recrind;

handles.Mtime=recMtime;
handles.Mtime2=recMtime2;
handles.Gmc=recG_mc;
handles.GMc=recG_Mc;

handles.S=recS;
handles.R=recR;
handles.D=recD;
handles.A=recA;
handles.J=recJ;

handles.Gstar=recGstar;
handles.F0=recF0;
handles.F1=recF1;
handles.F=recF;
handles.Psi=recPsi;


handles.nnS=nnS;
handles.nnR=nnR;
handles.nnD=nnD;
handles.nnA=nnA;
handles.nnJ=nnJ;

handles.NumDbyeachROI = recNumDcell;

    handles.mT=recmT;
    handles.MT=recMT;
    handles.numhL= rechLinc;
    
disp("M calculation end")

% disp("souka")
% nanmean(souka)
% sqrt(nanvar(souka))
% disp("soujou")
% nanmean(soujou)
% sqrt(nanvar(soujou))

datas=handles;
end
