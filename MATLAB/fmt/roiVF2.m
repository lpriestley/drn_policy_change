function out = roiVF2(Data, options, varargin)

% roiVF2 very similar to roiVF, but takes 'options' structure as input
%
% EDITED version of ViewFMRI, to allow roi selection
%
% Tool for viewing 4D FMRI datasets
%
% out = roiVF2(Data)
%
% out = roiVF2(Data,options,varargin)
%
% options.mean determines the data shown in the slice-by-slice lightbox.
% options.mean='m' shows the mean volume, this is the default.
% Can alternatively be a matrix volume of the same size as the data
%
% options.func is function that runs on middle mouse click
% default is roiplotseries, which just plots the time series
%
% options.funcd is function that runs when d is pressed
% default is side, which just plots the time series, with the ROI map
%
% varargin are any arguments to be passed to func and funcd
% The default 'func' (i.e. roiplotseries) takes time values as it's extra
% input, allowing non-uniformly spaced data to be plotted.

if nargin<2 options = []; end;

if isfield(options,'mean')  Mean=options.mean; else Mean = 'm'; end
if isfield(options,'func')  func=options.func; else func = 'roiplotseries'; end
if isfield(options,'funcd') funcd=options.funcd; else funcd = 'side'; end
if isfield(options,'cLims') cLims=options.cLims; else cLims = [min(min3(Data)) max(max3(Data))]; end
if isfield(options,'mask') mask=options.mask; else mask = []; end

% Rotate data to make imaging easier:
% Data = permute(Data,[2 1 3 4]);

warning off
if size(Data,1) ~= size(Data,2)
    fprintf('Image matrix has been padded with zeros to make it square');
    d1 = max([size(Data,1),size(Data,2)]);
    Data2 = zeros(d1,d1,size(Data,3),size(Data,4));
    Data2(1:size(Data,1),1:size(Data,2),:,:) = Data;
    Data = Data2;
    clear Data2;
end;

res = size(Data,1);

% Data = Data(:,end:-1:1,:,:);
% Data = swapmr_out(Data);

Index=0;
doFit = 0;

if (nargin<2),
    Mean='m';
end;

if length(size(Data))==2
    NumSlices=1;
end;


if ~isstr(Mean),
%     MeanImg=swapmr_out(Mean);
    MeanImg = Mean;

    if size(MeanImg,1) ~= size(MeanImg,2)
        d1 = max([size(MeanImg,1),size(MeanImg,2)]);
        MeanImg2 = zeros(d1,d1,size(MeanImg,3));
        MeanImg2(1:size(MeanImg,1),1:size(MeanImg,2),:) = MeanImg;
        MeanImg = MeanImg2;
        clear MeanImg2;
    end;

else,
    MeanImg=mean(Data,4);
end;


NumSlices=size(MeanImg,3);

if NumSlices>21
    AllData = MeanImg;
    MeanImg = AllData(:,:,Index+1:Index+21);
end
MeanImg=reshape(MeanImg,1,prod(size(MeanImg)));

if NumSlices<21
    tmp=zeros(1,size(Data,1)*size(Data,2)*21);
    tmp(1:size(Data,1)*size(Data,2)*NumSlices)=MeanImg;
    MeanImg=tmp;
end;

MeanImg=reshape(MeanImg,res,res,21);

LeftPanel=zeros(7*res,3*res);
for ct1=0:2
    for ct2=0:6
        LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)= flipud(MeanImg(:,:,1+ct1+ct2*3)');
    end;
end;

HndlImg=figure;
set(HndlImg,'Position',[    59   139   912   791]);

factor=2;

Xrange=[1,res];
Yrange=[1,res];
scale=res;
ImgNum=1;Xind=1;Yind=1;

LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
imagesc(LeftPanel);
axis('equal')
axis('off')
if NumSlices > 21
    title(sprintf('n for next 21 slices\np for previous 21 slices'));
end;

tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
MeanInt=MeanImg(Xind,Yind,ImgNum+Index);
subplot('position',[0.37,0.05,0.61,0.29]);
% feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});

RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
imagesc((MeanImg(:,:,ImgNum))',cLims)
set(RightImg,'XLim',Xrange,'YLim',Yrange,'YDir','normal');
axis('equal')
axis('tight')
set(RightImg,'XColor','red','YColor','red');
colormap(gray)

%--- ROI bit
maskNo = 1;
maskColor = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; .7 0 1; 1 .6 0; 0 .6 1]; %; 1 1 0];
if isempty(mask)
    noMaskColors = size(maskColor,1);
    mask = zeros(size(Data,1), size(Data,2), size(Data,3), noMaskColors);
else
    noMaskColors = size(mask,4);
    redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)
    if (noMaskColors<8), mask(:,:,:,8) = 0; end
end
dummyMask = zeros(size(mask)); zeroMask = zeros(size(mask));
%---

while 1==1

    TextImg=subplot('position',[0.8,0.41,0.18,0.55]);
    plot([1])
    cla
    axis('off')
    hold on
    set(TextImg,'XLim',[0,2],'YLim',[0,2]);
    text(0.3,1.8,sprintf('Slice No: %2d',ImgNum+Index));
    text(0.3,1.7,sprintf('       X: %2d  Y: %2d',Yind,Xind));
    text(0.3,1.6,sprintf('Intensity: %5.3f',MeanInt));
    line([0.05,1],[1.55,1.55],'Color','black')
    text(0.05,1.5,sprintf('+,-: zoom'));
    text(0.05,1.4,sprintf('<left>,<right>: select, deselect'));
    text(0.05,1.3,sprintf('<centre>: check voxel intensity '));
    text(0.05,1.2,sprintf('q,Q: quit'));
    text(0.05,1.1,sprintf('d,D: run ''funcd'' in new window'));
    text(0.05,1.0,sprintf('t,T: toggle masking'));
    text(0.05,0.9,sprintf('c,C: clear current mask'));
    text(0.05,0.8,sprintf('a,A: clear all masks'));
    text(0.05,0.7,sprintf('r,R/y,Y: +/- mask no.'));
    text(0.05,0.5,sprintf('Current mask: %d',maskNo));
    plot(1.5,.5,'s','MarkerSize',20,'Color',maskColor(maskNo,:),'MarkerFaceColor',maskColor(maskNo,:))
    
    plot(1,.15,'s','MarkerSize',100,'Color','white','MarkerFaceColor','white')
    
    posX = [.7 .7 .7 .7 1.3 1.3 1.3 1.3];
    posY = [.3 .2 .1 0 .3 .2 .1 0];
    for iColor = 1:noMaskColors
        text(posX(iColor),posY(iColor),num2str(sum3(mask(:,:,:,iColor))),'Color',maskColor(iColor,:))
    end
    
    drawflag=1;

    subplot(RightImg);
    [Y_c,X_c,B_c]=ginput(1);

    if gca==LeftImg

        dummy=Y_c; Y_c = X_c; X_c=dummy;
        
        tmp=1+floor(X_c/res)+3*floor(Y_c/res);
        if tmp<=NumSlices-Index
            ImgNum=tmp;
        else
            drawflag=0;
        end;

        if drawflag
            tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));

            MeanInt=MeanImg(Xind,Yind,ImgNum);
            subplot('position',[0.37,0.05,0.61,0.29]);
            feval(func, Data, mask, maskColor, doFit, varargin{:});
         
            redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)

            LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
            imagesc(LeftPanel);
            if NumSlices > 21
                title(sprintf('n for next 21 slices\np for previous 21 slices'));
            end;
            axis('off')
            axis('equal')
            Corner1=1+res*floor(X_c/res);
            Corner2=1+res*floor(Y_c/res);
            line([Corner1,Corner1],[Corner2,Corner2+res-1]);
            line([Corner1+res-1,Corner1+res-1],[Corner2,Corner2+res-1]);
            line([Corner1,Corner1+res-1],[Corner2,Corner2]);
            line([Corner1,Corner1+res-1],[Corner2+res-1,Corner2+res-1]);
        end;
%         title(['X = ' num2str(X_c) ', Y = ' num2str(Y_c)])

    end;

    if gca==RightImg
        if (B_c==2)|(B_c==double('+'))|(B_c==double('-'))
            if (round(X_c)>0)&(round(X_c)<=res)
                Xind=round(X_c);
            end;
            if (round(Y_c)>0)&(round(Y_c)<=res)
                Yind=round(Y_c);
            end;
        end;
        if B_c==2
            TSImg=subplot('position',[0.37,0.05,0.61,0.29]);
            cla
            tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
            MeanInt=MeanImg(Xind,Yind,ImgNum);
            feval(func, Data, mask, maskColor, doFit, varargin{:})
%             hold on
%             plot(tseries,'--');
        end;
        if (B_c==double('+'))&(scale>3)
            scale=round(scale/factor);
            Xrange=[Yind-0.5*scale,Yind+0.5*scale-1];
            Yrange=[Xind-0.5*scale,Xind+0.5*scale-1];
            %%Xrange=Xrange+sum(Xrange<1)-sum(Xrange>res);
            %%Yrange=Yrange+sum(Yrange<1)-sum(Yrange>res);
            Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
            Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);
            %%fprintf('X:%d Y:%d Xr: %d-%d Yr: %d-%d\n',Xind,Yind,min(Xrange),max(Xrange),min(Yrange),max(Yrange))
            
            redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)

        end;
        if (B_c==double('-'))&(scale < res-1)
            scale=round(scale*factor);
            Xrange=[Yind-0.5*scale,Yind+0.5*scale-1];
            Yrange=[Xind-0.5*scale,Xind+0.5*scale-1];
            Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
            Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);

            redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)

        end;
        
        if (B_c==double('y'))|(B_c==double('Y'))
            maskNo = maskNo + 1;
            if maskNo > noMaskColors
                maskNo = 1;
            end
        end
        if (B_c==double('r'))|(B_c==double('R'))
            maskNo = maskNo - 1;
            if maskNo < 1
                maskNo = noMaskColors;
            end
        end
        
        if (B_c==1)
            if (round(X_c)>0)&(round(X_c)<=res)
                Xind=round(X_c);
            end;
            if (round(Y_c)>0)&(round(Y_c)<=res)
                Yind=round(Y_c);
            end;
            if ~mask(Yind,Xind,ImgNum+Index,find(1:noMaskColors ~= maskNo))
                % only alter the mask if another mask hasn't already used that
                % voxel
                if mask(Yind,Xind,ImgNum+Index,maskNo) == 0
                    mask(Yind,Xind,ImgNum+Index,maskNo) = 1;
                end
            end
            redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)
            
            subplot('position',[0.37,0.05,0.61,0.29]);
            cla
            feval(func, Data, mask, maskColor, doFit, varargin{:});
        end

        if (B_c==3)
            if (round(X_c)>0)&(round(X_c)<=res)
                Xind=round(X_c);
            end;
            if (round(Y_c)>0)&(round(Y_c)<=res)
                Yind=round(Y_c);
            end;
            if ~mask(Yind,Xind,ImgNum+Index,find(1:noMaskColors ~= maskNo))
                % only alter the mask if another mask hasn't already used that
                % voxel
                if mask(Yind,Xind,ImgNum+Index,maskNo) == 1
                    mask(Yind,Xind,ImgNum+Index,maskNo) = 0;
                end
            end
            redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)

            subplot('position',[0.37,0.05,0.61,0.29]);
            cla
            feval(func, Data, mask, maskColor, doFit, varargin{:});
        end
        
    end;
    
    if (B_c==double('T'))|(B_c==double('t'))
        dummyMask = mask;
        mask = zeroMask;
        zeroMask = dummyMask;
        redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)
    end
    
    
    if (B_c==double('q'))|(B_c==double('Q')) % i.e. quit!
%         out.mask = permute(mask,[2 1 3 4]);
        out.mask=mask;
        out.maskColor = maskColor;
        break
    end;
    if (B_c==double('d'))|(B_c==double('D'))
        figure
        if strcmp(funcd,'side')

%             set(gcf,'Position',[            49         329        1195         421])
            set(gcf,'Position',[ 104   549   879   250])
            LeftFig = subplot(121);
            RightFig = subplot(122);
            set(LeftFig, 'Position',[0.0134    0.0880    0.2346    0.8490]);
            set(RightFig,'Position',[0.2960    0.1640    0.6628    0.7120]);

            subplot(LeftFig)
            DisplayImg = Intense2rgb(MeanImg(:,:,ImgNum)');
            for j = 1:noMaskColors
                [y,x] = find(mask(:,:,ImgNum+Index,j));
                for i = 1:length(x)
                    DisplayImg(x(i),y(i),:) = maskColor(j,:);
                end
            end
            image(DisplayImg)
            axis equal tight
            set(gca,'xtick',[],'ytick',[],'YDir','normal')

%             title(['Slice ',num2str(ImgNum+Index)]);

            subplot(RightFig)
            feval(func, Data, mask, maskColor, doFit, varargin{:})
            allText   = findall(gcf, 'type', 'text');
            allAxes   = findall(gcf, 'type', 'axes');
            allFont   = [allText; allAxes];
            set(allFont,'FontSize',12,'Fontweight','bold')
        else
            feval(funcd, Data, mask, maskColor, doFit, varargin{:});
        end
        figure(HndlImg);
    end;
    if (B_c==112)|(B_c==80)
        %P or p
        if (Index>0)
            Index=Index-21;
            MeanImg = zeros(res,res,21);
            MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices]));
            LeftPanel=zeros(7*res,3*res);
            for ct1=0:2
                for ct2=0:6
                    LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=flipud(MeanImg(:,:,1+ct1+ct2*3)');
                end;
            end;
            subplot(LeftImg);
            imagesc(LeftPanel);
            axis('equal')
            axis('off')
            if NumSlices > 21
                title(sprintf('n for next 21 slices\np for previous 21 slices'));
            end;
        end;
    end;
    if (B_c==110)|(B_c==78)
        %N or n
        if ((Index+21)<NumSlices)
            Index=Index+21;
            MeanImg = zeros(res,res,21);
            MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices]));
            LeftPanel=zeros(7*res,3*res);
            for ct1=0:2
                for ct2=0:6
                    LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=flipud(MeanImg(:,:,1+ct1+ct2*3)');
                end;
            end;
            subplot(LeftImg);
            imagesc(LeftPanel);
            axis('equal')
            axis('off')
            if NumSlices > 21
                title(sprintf('n for next 21 slices\np for previous 21 slices'));
            end;
        end;
    end;
    
    if (B_c==double('f'))|(B_c==double('F'))
        if doFit  doFit = 0; else  doFit = 1;  end
    end
    
    if (B_c==double('a'))|(B_c==double('A'))
        mask = zeros(size(mask));
        redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims);
    end

    if (B_c==double('c'))|(B_c==double('C'))
       mask(:,:,:,maskNo) = zeros(res,res,size(Data,3)); 
       redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims);
    end
    
end;

%%%%%%%%%%%%%%%%%%%%%%
function redrawRightImg(MeanImg,ImgNum,Index,Xrange,Yrange,mask,maskColor,noMaskColors,maskNo,cLims)
    RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
    DisplayImg = Intense2rgb(MeanImg(:,:,ImgNum)',[],cLims);
    for j = 1:noMaskColors
        [y,x] = find(mask(:,:,ImgNum+Index,j));
        for i = 1:length(x)
            DisplayImg(x(i),y(i),:) = maskColor(j,:);
        end
    end
    image(DisplayImg)
    axis equal
    set(RightImg,'XLim',Xrange,'YLim',Yrange,'YDir','normal');
    title(['Slice ',num2str(ImgNum+Index)]);
    set(RightImg,'XColor','red','YColor','red');

%%%%%%%%%%%%%%%%%%%%  
function roiplotseries(Data, mask, maskColor, doFit, TIs)
if nargin < 5
    TIs = 1:size(Data,4);
end

hold on
for i = 1:size(maskColor,1)
    noPixInMask = sum3(mask(:,:,:,i));
    
    if noPixInMask == 1
        yData = mean(mask_series3(Data,mask(:,:,:,i)),1);
        plot(TIs,yData,'-o','Color',maskColor(i,:),'LineWidth',2)
    elseif noPixInMask > 1
        yData = mean(mask_series3(Data,mask(:,:,:,i)),1);
        yStdE = std(mask_series3(Data,mask(:,:,:,i)),0,1); yStdE = yStdE/sqrt(length(yStdE));
        errorbar(TIs,yData,yStdE,'-o','Color',maskColor(i,:),'LineWidth',2)
    end
end
grid on
hold off;

%%%%%%%%%%%%%%%%%%%%  
function roiplotseries_noline(Data, mask, maskColor, doFit, TIs)
if nargin < 5
    TIs = 1:size(Data,4);
end

hold on
for i = 1:size(maskColor,1)
    noPixInMask = sum3(mask(:,:,:,i));
    
    if noPixInMask == 1
        yData = mean(mask_series3(Data,mask(:,:,:,i)),1);
        plot(TIs,yData,'x','Color',maskColor(i,:),'LineWidth',2)
    elseif noPixInMask > 1
        yData = mean(mask_series3(Data,mask(:,:,:,i)),1);
        yStdE = std(mask_series3(Data,mask(:,:,:,i)),0,1); yStdE = yStdE/sqrt(length(yStdE));
        errorbar(TIs,yData,yStdE,'x','Color',maskColor(i,:),'LineWidth',2)
    end
end
grid on
hold off;
