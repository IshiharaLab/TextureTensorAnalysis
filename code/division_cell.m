function data=division_cell(handles)

NumberImages = handles.In;

numD = cell(NumberImages,1);
numA = cell(NumberImages,1);
Divtime = cell(NumberImages-1,1);
numD{1}= [];
BW=handles.AreaLabel{1};
rec_out_cell = handles.outcell;
[Lm,Ln]=size(BW);
for inum=2:NumberImages
    xa=handles.pos{inum}(:,1);
    ind2=handles.ind2{inum-1};
    out_cell = rec_out_cell{inum};

    tmpL=1:length(xa);
    ca=tmpL(~ismember(tmpL,ind2));  % before connection nothing
    in_center = tmpL(~ismember(tmpL,out_cell)); % eject most outer cell in tissue
    numD{inum}= in_center(ismember(in_center,ca));

    % apothosis
    ind1=handles.ind1{inum-1};
    xb=handles.pos{inum-1}(:,1);
    tmpL=1:length(xb);
    ca=tmpL(~ismember(tmpL,ind1));  % before connection nothing
    out_cell = rec_out_cell{inum-1};
    in_center = tmpL(~ismember(tmpL,out_cell)); % eject most outer cell in tissue
    numA{inum-1}= in_center(ismember(in_center,ca));
end
numA{NumberImages}= [];
handles.numD = numD;    % nothing to connect with the previous frame.
handles.numA = numA;    % nothing to connect with the previous frame.

D_num=cell(5,NumberImages);
for i=1:5
    D_num{i,1}=[];
end

for inum=1:NumberImages-1
    BArea=handles.AreaLabel{inum};
    Area=handles.AreaLabel{inum+1};
    BW = handles.I(:,:,inum+1);

    r0=handles.pos{inum+1};
    nn=numD{inum+1};
    cnnt=handles.nearC{inum+1};

    ind1=handles.ind1{inum};
    ind2=handles.ind2{inum};

    if isempty(nn)    % if there arenot division on inum+1, D_num{inum+1} can be decided.
        for i=1:5
            if isempty(D_num{i,inum})
                D_num{i,inum+1} = [];
            else
                D_num{i,inum+1} = ind2(ismember(ind1,D_num{i,inum}));
            end
        end
    else
        % D de connect-to-original (nn2) and new apear (nn)
        nn2=[];
        nnRec =nn;
        for i=1: length(nn)
            link=logical(sum(cnnt==nn(i),2));
            numLink=sum(link);
            A=cnnt(link,:);
            nC=A(A~=nn(i));

            tmp=[];
            for j=1:numLink
                tmp=[tmp  sum(ismember(BArea,ind1(ind2==nC(j))),'all') - sum(ismember(Area,nC(j)),'all')];
            end
            if max(tmp)>0
                if sum(tmp==max(tmp))~=1
                    tmC=1:length(tmp);
                    tma = tmC(tmp==max(tmp));
                    nn2 = [nn2 nC(tma(1))];
                else
                    nn2 = [nn2 nC(tmp==max(tmp))];
                end
            else    %  all  cells are increased area.
                nnRec=nnRec(nnRec~=nn(i));
            end
        end
        nn =nnRec;

        Lnn2 = length(nn2);
        nn0=[];
        for iii=1:Lnn2
            nn0 = [nn0; ind1(ismember(ind2,nn2(iii)))];
        end


        Divtime{inum}=[nn0,nn',nn2'];
        for i=1:5
            DDD=[];
            if ~isempty(D_num{i,inum})
                DDD=D_num{i,inum};
                for jj=1:length(nn2)
                    tmp=ind1(ind2==nn2(jj));
                    DDD = DDD(~ismember(DDD,tmp));
                end
                DDD = ind2(ismember(ind1,DDD));
            end
            D_num{i,inum+1} =DDD';
        end
        % At this point, the divided cell that has already been counted was selected.
        for jj=1:length(nn2)
            tmp=ind1(ind2==nn2(jj));
            idx=0;
            %?��?Examining whether new division cells are covered with D_num.
            for ii=1:4
                if idx==0
                    if sum(ismember(D_num{ii,inum},tmp))~=0
                        D_num{ii+1,inum+1}=[D_num{ii+1,inum+1} nn(jj) nn2(jj)];
                        idx=1;
                    end
                end
            end
            if idx==0
                if sum(ismember(D_num{5,inum},tmp))~=0
                    D_num{5,inum+1}=[D_num{5,inum+1} nn(jj) nn2(jj)];
                else
                    D_num{1,inum+1}=[D_num{1,inum+1} nn(jj) nn2(jj)];
                end
            end
        end
    end
end


disp("division end")

handles.numDTot = D_num;
handles.Divtime = Divtime;

data = handles;

end