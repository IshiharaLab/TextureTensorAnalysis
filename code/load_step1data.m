function data=load_step1data(str)
    load(str,'data')
    
    BW = data.BWI;
    NumberImages = size(BW,3);
    data.In = NumberImages;
    data.in = 1; % current image number
    data.I = BW;
end