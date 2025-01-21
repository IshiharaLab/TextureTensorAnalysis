function [M2,num2] = Mtensor(M,num,sel,cnnt,xq,yq,vt4)
    link = logical(sum(cnnt==sel,2));
    tmpM = cnnt(link,:);
    each_ind = tmpM(tmpM~=sel);
    tmp=length(each_ind);
    lx = xq(each_ind)-xq(sel);
    ly = yq(each_ind)-yq(sel);
    M2 = 0.5*[sum(lx.*lx), sum(lx.*ly); sum(lx.*ly), sum(ly.*ly)];
    num2 = num + length(lx);
    % M+0.5*[ ] % added 4vtx
    chk4 = sum(ismember(vt4,sel),1:2);
    if chk4 ~= 0
        a = logical(sum(ismember(vt4,sel),2));
        tmp4vtx = vt4(a,:);
        each4vt = [each_ind;sel];
        tmp2 = tmp4vtx(~ismember(tmp4vtx,each4vt));
        tmp=tmp+0.5*length(tmp2);
        lx = xq(tmp2)-xq(sel);
        ly = yq(tmp2)-yq(sel);
        M2 = M2 + 0.25*[sum(lx.*lx), sum(lx.*ly); sum(lx.*ly), sum(ly.*ly)];
        num2 = num + 0.5*length(lx);
    end
    M2 = M + M2/tmp;
end