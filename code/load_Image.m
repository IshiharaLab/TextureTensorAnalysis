function data=load_Image(inptdata,name1,name2)
    
    data=inptdata;
    
    chkk = imread(name1);
    if isa(chkk,'uint8')~=1
        errordlg('Please select the 8bit type')
        return
    end
        
    InfoImage = imfinfo(name1);
    mImage = InfoImage(1).Width;
    nImage = InfoImage(1).Height;
    NumberImages = length(InfoImage);
    Image = zeros(nImage,mImage,NumberImages); % image(red)
    BWImage = zeros(nImage,mImage,NumberImages); % BW channel
    
    h=msgbox('Now reading ...','modal');
    for i=1:NumberImages
        temp = imread(name1,'Index',i,'Info',InfoImage);
        Image(:,:,i) = temp(:,:,1);
        BB=Image(:,:,i)/255;
        BW = imbinarize(BB(:,:));
        BWImage(:,:,i) = BW(:,:);
    end
    
    close(h);
    
    
    tic
    
    [m,n]=size(BWImage(:,:,1));
    BWImage(1,:,:)=0;
    BWImage(:,1,:)=0;
    BWImage(m,:,:)=0;
    BWImage(:,n,:)=0;
    BWImage = 1-BWImage;
    
    
    disp('=============')
    disp('1.5 Load PIV data')

    % load the PIV-data
     load(name2);

     xpiv=x;
     ypiv=y;

    data.pivx = xpiv;
    data.pivy = ypiv;
    data.pivu = u_filtered;
    data.pivv = v_filtered;

    se = strel('diamond',2);
    x_rec=cell(NumberImages,1);
    nearcell_rec=cell(NumberImages,1);
    num_cell_mat=cell(NumberImages,1);
    vertex4_rec=cell(NumberImages,1);
    Outregion = cell(1,NumberImages);
    
    disp('=============')
    disp('1.7 neighboring cell')

    parfor j=1:NumberImages
        BW=BWImage(:,:,j);
        L = bwlabeln(BW,4);
        Mx=max(L(:));
    
        Stmp=0.0;
        Sind=0; 
        for k=1:Mx
            if( Stmp < sum(L(:)==k) )
                Stmp=sum(L(:)==k);
                Sind=k; 
            end         
        end
        Outregion{j} = (L==Sind);
        L(L(:,:)==Sind)=0;
        L(L(:,:)>Sind) = L(L(:,:)>Sind)-1;
        num_cell_mat{j}=L;
        
        x=[];
        for i=1:Mx
            [tmpx,tmpy]=find(L==i);
            if ~isnan(tmpx)
                x=[ x; mean(tmpx), mean(tmpy) ];
            end
        end
       
        [m,n]=size(L);
        [a1,b1]=find(L(:,:)==0);
        b2=b1(a1~=1 & b1~=1 & a1~=m & b1~=n);
        a2=a1(a1~=1 & b1~=1 & a1~=m & b1~=n);
        ll=length(a2);
        nearcell=[];
        vertex4=[];
        for j3=1:ll
            tmp=unique(L(a2(j3)-1:a2(j3)+1,b2(j3)-1:b2(j3)+1))';
            if length(tmp(tmp~=0))==2
                nearcell=[nearcell;tmp(tmp~=0)];
            elseif length(tmp(tmp~=0))==4
                vertex4=[vertex4; tmp(tmp~=0)];
            end
        end
        nearcell2 = unique(nearcell,'rows');
         
         x_rec{j}=x;
         nearcell_rec{j}=nearcell2;
         vertex4_rec{j}=vertex4;
    end
    
        disp('=============')
        disp('1.8 Outer cell')
    % Outside cell
    rec_out_cell = cell(1,NumberImages);
    parfor inum = 1:NumberImages
        M2 = Outregion{inum};
        [lm,ln] = size(M2);
        Mrec = zeros(lm,ln);
        M = zeros(lm,ln);
        tmp = M2(2:lm,1:ln)-M2(1:lm-1,1:ln);
        Mrec(1:lm-1,1:ln) = (tmp==1);
        M(2:lm,1:ln) = (tmp==-1);
        Mrec = M+Mrec;
        
        tmp = M2(1:lm,2:ln)-M2(1:lm,1:ln-1);
        M(1:lm,1:ln-1) = (tmp==1);
        Mrec = M+Mrec;
        M(1:lm,2:ln) = (tmp==-1);
        Mrec = M+Mrec;
        M = (Mrec~=0);
        
        L = num_cell_mat{inum};
        [m,n] = size(L);
        [a2,b2] = find(M(:,:)==1);
        
        ll=length(a2);
        edgecell=[];
        for j3=1:ll
            tmp=L(a2(j3)-1:a2(j3)+1,b2(j3)-1:b2(j3)+1);
            if length(unique(tmp(~ismember(tmp,0))))==1
                edgecell = [edgecell unique(tmp(~ismember(tmp,0)))];
            end
        end
        rec_out_cell{inum}=unique(edgecell);
        
    end

    %%
    % ffnam = strcat('pre_pos_link_step1.mat');
    % filename=ffnam;
    
    data.pos=x_rec;
    data.BWI=BWImage;
    data.nearC=nearcell_rec;
    data.AreaLabel=num_cell_mat;
    % data.pivx = xpiv;
    % data.pivy = ypiv;
    % data.pivu = u_filtered;
    % data.pivv = v_filtered;
    data.outcell = rec_out_cell;
    data.Outregion=Outregion;
    data.vtx4 = vertex4_rec;
    
    % h=msgbox('Now saving ...','modal');
    % m = matfile(filename);
    % save(filename,'data')
    % close(h);


    rmax = 400.0;
    
    rec_ind1 = cell(NumberImages-1,1);
    rec_ind2 = cell(NumberImages-1,1);
    rec_rpiv = cell(NumberImages-1,1);
    
    tmp_r1 = cell(NumberImages-1,1);
    tmp_r2 = cell(NumberImages-1,1);
    for inum=1:NumberImages-1
        tmp_r1{inum} =  x_rec{inum};
        tmp_r2{inum} =  x_rec{inum+1};
    end
    
    disp('=============')
    disp('1.9 Tracking')


    BW=BWImage(:,:,1);
    parfor inum=1:NumberImages-1
        % %PIV collection
        r1 = tmp_r1{inum};
        R2 = tmp_r2{inum};
        k=convhull(r1(:,1),r1(:,2));
        [in,~] = inpolygon(xpiv{inum},ypiv{inum},r1(k,2),r1(k,1));
        
        vf = data.pivv{inum};
        uf = data.pivu{inum};
        R1=zeros(length(r1),2);
        for num_cell=1:length(r1)
            locPIV=(xpiv{inum}>(r1(num_cell,2)-20) & xpiv{inum}<(r1(num_cell,2)+20) & ypiv{inum}>(r1(num_cell,1)-20) & ypiv{inum}<(r1(num_cell,1)+20));
            locPIV2=in+locPIV==2;
            if sum(locPIV2,1:2) ==0
                R1(num_cell,:) = r1(num_cell,:);
            else
                 R1(num_cell,:) = r1(num_cell,:)+mean([vf(locPIV2),uf(locPIV2)]);
            end
        end
        rec_rpiv{inum}=R1;
        
        D = sqdist(R1',R2');
        Dm = 1.05*max(D(:));
        D(D(:) >= rmax) = Inf;
        [assignments, ~,~] = assignjv(D,140);
        rec_ind1{inum}=assignments(:,1);
        rec_ind2{inum}=assignments(:,2);
        
    end
    
    toc

    
    ffnam = strcat('pos_link_step1.mat');
    [filename, PathName]  = uiputfile(ffnam,'Save file name (.mat)');
    if filename == 0
        return
    end
    
    data.ind1=rec_ind1;
    data.ind2=rec_ind2;
    data.rpiv=rec_rpiv;
    % data.pos=x_rec;
    % data.BWI=BWImage;
    % data.nearC=nearcell_rec;
    % data.AreaLabel=num_cell_mat;
    % data.pivx = xpiv;
    % data.pivy = ypiv;
    % data.pivu = u_filtered;
    % data.pivv = v_filtered;
    % data.outcell = rec_out_cell;
    % data.vtx4 = vertex4_rec;
    % data.Outregion=Outregion;
    
    
    h=msgbox('Now saving ...','modal');
    m = matfile(filename);
    save(filename,'data')
    close(h);


end
