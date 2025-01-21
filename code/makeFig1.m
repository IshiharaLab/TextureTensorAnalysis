function makeFig1(handles)

NumberImages=handles.In;

figh=figure;
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile

inum = 1;
roi      = handles.roi{inum};
d_roi        = size(roi,2);
mask_sepa=handles.IROI(:,:,inum);

BW = handles.I(:,:,inum);
[LLm,LLn]=size(BW);

position=[0.95*LLn 0.1*LLm];
value=[inum];
BW2 = insertText(BW,position,value,'AnchorPoint','LeftBottom');

L=handles.AreaLabel{inum};
OCell=handles.outcell{inum};
mask_sepa=handles.IROI(:,:,inum);
RGB2 = label2rgb(mask_sepa,'pink','c','shuffle');
[lm,ln,col]=size(RGB2);
clear mask_sepa;
OCI=ismember(L,OCell);
I3 = uint8(zeros(lm,ln,col));
A=RGB2(:,:,1);
B=RGB2(:,:,2);
C=RGB2(:,:,3);
A(OCI==1)=210;
B(OCI==1)=210;
C(OCI==1)=210;
I3(:,:,1)=A;
I3(:,:,2)=B;
I3(:,:,3)=C;
imshow(I3)
hold on

recMtime= handles.Mtime;
Mrec=recMtime{inum};

xcent = handles.roipos{inum};
tmpt=0:10:360;
for jroi=1:d_roi
    [Vec,Dval] = eig(Mrec{jroi});
    Dval=0.2*Dval;
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
title("first frame")


nexttile


inum = NumberImages;
roi      = handles.roi{inum};
d_roi        = size(roi,2);
mask_sepa=handles.IROI(:,:,inum);

BW = handles.I(:,:,inum);
[LLm,LLn]=size(BW);

position=[0.95*LLn 0.1*LLm];
value=[inum];
BW2 = insertText(BW,position,value,'AnchorPoint','LeftBottom');

L=handles.AreaLabel{inum};
OCell=handles.outcell{inum};
mask_sepa=handles.IROI(:,:,inum);
RGB2 = label2rgb(mask_sepa,'pink','c','shuffle');
[lm,ln,col]=size(RGB2);
clear mask_sepa;
OCI=ismember(L,OCell);
I3 = uint8(zeros(lm,ln,col));
A=RGB2(:,:,1);
B=RGB2(:,:,2);
C=RGB2(:,:,3);
A(OCI==1)=210;
B(OCI==1)=210;
C(OCI==1)=210;
I3(:,:,1)=A;
I3(:,:,2)=B;
I3(:,:,3)=C;
imshow(I3)
hold on

recMtime= handles.Mtime;
Mrec=recMtime{inum};

xcent = handles.roipos{inum};
tmpt=0:10:360;
for jroi=1:d_roi
    [Vec,Dval] = eig(Mrec{jroi});
    Dval=0.2*Dval;
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
title("last frame")

str=['Figure/exampleFig_M_and_ROI.jpg'];
saveas(figh,str)
close(figh);
end