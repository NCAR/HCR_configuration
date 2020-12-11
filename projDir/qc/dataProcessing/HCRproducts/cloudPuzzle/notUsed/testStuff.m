%% Two thresholds, find inner and outer boundaries

BW(reflMap>=thresh12(1) & reflMap<thresh12(2))=1;
        
        [B,L,N,A] = bwboundaries(BW);
        
        cutBound=100;
        cutBoundInner=1000;
        goodB={};
        figure; imshow(BW); hold on;
        % Loop through object boundaries
        for k = 1:N
            % Boundary k is the parent of a hole if the k-th column
            % of the adjacency matrix A contains a non-zero element
            if (nnz(A(:,k)) > 0)
                %                 boundary = B{k};
                %                 plot(boundary(:,2),...
                %                     boundary(:,1),'r','LineWidth',2);
                % Loop through the children of boundary k
                keepParent=1;
                for l = find(A(:,k))'
                    boundary = B{l};
                    if size(boundary,1)>cutBoundInner
                        goodB{end+1}=boundary;
                        keepParent=0;
                        
                        plot(boundary(:,2),...
                            boundary(:,1),'g','LineWidth',2);
                    end
                end
                if keepParent==1
                    boundary = B{k};
                    if size(boundary,1)>cutBound
                        plot(boundary(:,2),...
                            boundary(:,1),'r','LineWidth',2);
                        goodB{end+1}=boundary;
                    end
                end
            else
                boundary = B{k};
                if size(boundary,1)>cutBound
                    plot(boundary(:,2),...
                        boundary(:,1),'r','LineWidth',2);
                    goodB{end+1}=boundary;
                end
            end
        end
        
        %% Adaptive thresholding
        
               
        
                % Adjust data to span data range.
                Iadj = imadjust(I);
        
                % Threshold image - adaptive threshold
                BW = imbinarize(Iadj, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'dark');
                BW(isnan(reflMap))=0;
        
                BW = imfill(BW,'holes');
                BW = bwareaopen(BW,1000);
        
                threshMask=double(BW);
                threshMask(threshMask==0)=nan;
        
                properties = regionprops(BW, {'Area', 'Eccentricity', 'EquivDiameter',...
                    'Extent', 'FilledArea', 'MajorAxisLength', 'MinorAxisLength',...
                    'Orientation', 'Perimeter', 'Solidity'});
        
                imageRegionAnalyzer(BW);
        I = mat2gray(reflMap);
        I(reflMap==1)=0;
        I=imadjust(I);
        J = im2uint8(I);
        
        [L,Centers] = imsegkmeans(J,2);
        
        regionsOut=im2double(L);
        regionsOut(isnan(reflMap))=nan;
        
        %% Image morphology tool
       
        imageMorphology(I)
        
        %% K means clustering
        
        I = mat2gray(reflMap);
        I(reflMap==1)=0;
        I=imadjust(I);
        J = im2uint8(I);
        
        [L,Centers] = imsegkmeans(J,2);
        
        regionsOut=im2double(L);
        regionsOut(isnan(reflMap))=nan;
        
        function allPointsRefls= attachPoints(maskLabel,reflPadded,threshIn)
% maskLabel is a matrix created from a thresholded mask where each area has
% one number. We try to attach all remaining non nan reflectivity values to
% one of these areas by an iterative approach which also uses thresholding

allPointsRefls=nan(size(reflPadded));
reflSteps=fliplr(floor(min(min(reflPadded))):3:threshIn);

numAreas=max(max(maskLabel));
maskLabelNew=maskLabel;
close all
figure
surf(maskLabel,'edgecolor','none');
view(2)
colormap lines


for ii=2:length(reflSteps)
    oldMask=zeros(size(maskLabel));
    newMask=zeros(size(maskLabel));
    
    oldMask(~isnan(maskLabelNew))=1;    
    newMask(reflPadded>reflSteps(ii))=1;
    
    % Check which areas contain data from old mask and which are new
%     oldPixels=find(oldMask==1);
%         
%     newFilled=imfill(newMask,'holes');
%     newAreas=bwconncomp(newFilled);
%     newAreasPix=newAreas.PixelIdxList;
%     
%     for jj=1:newAreas.NumObjects
%         overlapping=intersect(oldPixels,newAreasPix{jj});
%         if isempty(overlapping)
%             numAreas=numAreas+1;
%             maskLabel(newAreasPix{jj})=numAreas;
%             oldMask(newAreasPix{jj})=1;
%             maskLabelNew(newAreasPix{jj})=numAreas;
%         end
%     end
        
    [oldR oldC]=find(oldMask==1);
    [addR addC]=find(oldMask==0 & newMask==1);
    %[rowA, colA] = find(oldMask==1);
idx = knnsearch([oldR oldC], [addR addC]);
%nearestR=oldR(idx);
%nearestC=oldC(idx);
nearest_OldValue = maskLabelNew(sub2ind(size(maskLabel), oldR(idx), oldC(idx)));
    maskLabelNew(sub2ind(size(maskLabel), addR, addC))=nearest_OldValue;
    
    figure
colormap lines

    surf(maskLabelNew,'edgecolor','none');
    view(2)
    
    
    
%     x=cat(1,oldR,oldC);
%     Mdl = KDTreeSearcher(x);
%     [addR addC]=find(oldMask==0 & newMask==1);
%     newpoint=cat(1,addR,addC);
%     [n,d] = knnsearch(Mdl,newpoint,'k',1);
%     
%     Idx = knnsearch(oldPixelsAgain,addPixels);
%     maskLabelNew(addPixels)=maskLabel(Idx);
    
end
    

    
end

    %% Smooth with convolution
    
    % Create mask for convolution
    radius=5;
    numPix=radius*2+1;
    [rr cc] = meshgrid(1:numPix);
    cirMask = sqrt((rr-(radius+1)).^2+(cc-(radius+1)).^2)<=radius;
    cirMask=double(cirMask);
    
    % Normalize
    cirMask=cirMask./(sum(reshape(cirMask,1,[])));
    
    % Convolution
    reflConv=nanconv(reflExt,cirMask);
    reflConv(isnan(reflExt))=nan;