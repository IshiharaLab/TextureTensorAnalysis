function datas=sepa_ROI(handles,sided)

%�@Minimum number of cells in ROI
cutROI = 7;


inum =1;

if isfield(handles,'xyROI')==1
    xy = handles.xyROI;
else
    areaA=handles.AreaLabel{inum};
    [Lm,Ln]=size(areaA);
    xy=[0 0;0 Lm;Ln Lm;Ln 0];
end
NumberImages = handles.In;
roiT=cell(NumberImages,1);
roiposT=cell(NumberImages,1);

Mxy=max(xy);
mxy=min(xy);

nd = round((Mxy(1)-mxy(1))/sided);
md = round((Mxy(2)-mxy(2))/sided);

sdx=(Mxy(1)-mxy(1))/nd;
sdy=(Mxy(2)-mxy(2))/md;
roi=cell(1,nd*md);
for j=1:nd
    for j2=1:md
        roi{md*(j-1)+j2}=[mxy(2)+sdy*(j2-1),mxy(1)+sdx*(j-1); mxy(2)+sdy*j2,mxy(1)+sdx*(j-1); mxy(2)+sdy*j2 ,mxy(1)+sdx*j; mxy(2)+sdy*(j2-1) , mxy(1)+sdx*j];
    end
end


idROI=cell(1,nd*md);
xcent=cell(1,nd*md);
xq       = handles.pos{inum}(:,2);
yq       = handles.pos{inum}(:,1);
oCell   = handles.outcell{inum};
numID=[];
d_roi=size(roi,2);
for j=1:d_roi
    [in,~] = inpolygon(xq,yq,xy(:,1),xy(:,2));
    c=1:size(in);
    inPoly = c(in);    
    [in,~] = inpolygon(xq,yq,roi{j}(:,2),roi{j}(:,1));
    c=1:size(in);
    inROI = c(in);    
    tmp = inPoly(ismember(inPoly,inROI));         % id in each ROI
    a = tmp(~ismember(tmp,oCell));         % eject the outer cell
    idROI{j} = a;         % eject the outer cell
    xcent{j} = [mean(xq(a)),mean(yq(a))];
    numID=[numID, length(a)];
end
aa=(numID>cutROI);
A=idROI(aa);
B=xcent(aa);
idROI=A;
xcent=B;
roiT{inum}=idROI;
roiposT{inum}=xcent;

d_roi= length(idROI);

areaA=handles.AreaLabel{inum};
mask_sepa=zeros(size(areaA));
for j=1:d_roi
    mask_sepa=mask_sepa+j*ismember(areaA,idROI{j});
end

% figure
% L=handles.AreaLabel{inum};
% OCell=handles.outcell{inum};
% RGB2 = label2rgb(mask_sepa,'pink','c','shuffle');
% [lm,ln,col]=size(RGB2);
% OCI=ismember(L,OCell);
% I3 = uint8(zeros(lm,ln,col));
% A=RGB2(:,:,1);
% B=RGB2(:,:,2);
% C=RGB2(:,:,3);
% A(OCI==1)=210;
% B(OCI==1)=210;
% C(OCI==1)=210;
% I3(:,:,1)=A;
% I3(:,:,2)=B;
% I3(:,:,3)=C;
% imshow(I3)

[Lm,Ln]=size(mask_sepa);
ImageROI = zeros(Lm,Ln,NumberImages);
ImageROI(:,:,inum) = mask_sepa; 
clear mask_sepa;

for inum = 1:NumberImages-1
    ind1    = handles.ind1{inum};
    ind2    = handles.ind2{inum};
    oCell   = handles.outcell{inum+1};
    xq       = handles.pos{inum+1}(:,2);
    yq       = handles.pos{inum+1}(:,1);

    if isempty(handles.Divtime{inum})
        cell_D=nan;
    else
        cell_D = handles.Divtime{inum};
        Db = handles.Divtime{inum}(:,1);
    end
    
    
    for j=1:d_roi
        tmp = idROI{j};
        link_next   = ind2(ismember(ind1,tmp));
        a1 = link_next(~ismember(link_next,oCell));         % eject the outer cell
        a12 = a1(:);
        if ~isnan(cell_D)
            Db_roi = tmp(ismember(tmp,Db));
            if ~isempty( ismember(Db,Db_roi))
                idD = ismember(Db,Db_roi);
                tmp2 = cell_D(idD,2:3);
                tmp3 = tmp2(:);
                a2 = tmp3(~ismember(tmp3,oCell));         % eject the outer cell
                a22 = a2(:);
            else
                a22=[];
            end
        else
            a22=[];
        end
        a3 = [a12;a22]; 
        a4 = unique(a3);
        idROI{j} = a4;          
        xcent{j} = [mean(xq(a4)),mean(yq(a4))];
    end
    roiT{inum+1}=idROI;
    roiposT{inum+1}=xcent;
    
areaA=handles.AreaLabel{inum+1};
mask_sepa=zeros(size(areaA));
for j=1:d_roi
    mask_sepa=mask_sepa+j*ismember(areaA,idROI{j});
end
ImageROI(:,:,inum+1) = mask_sepa; 
end

% figure
% L=handles.AreaLabel{inum};
% OCell=handles.outcell{inum};
% RGB2 = label2rgb(mask_sepa,'pink','c','shuffle');
% [lm,ln,col]=size(RGB2);
% OCI=ismember(L,OCell);
% I3 = uint8(zeros(lm,ln,col));
% A=RGB2(:,:,1);
% B=RGB2(:,:,2);
% C=RGB2(:,:,3);
% A(OCI==1)=210;
% B(OCI==1)=210;
% C(OCI==1)=210;
% I3(:,:,1)=A;
% I3(:,:,2)=B;
% I3(:,:,3)=C;
% imshow(I3)

handles.roi=roiT;
handles.roipos=roiposT;
handles.sided = sided;
handles.IROI = ImageROI; 

datas = handles;

disp("End")

end