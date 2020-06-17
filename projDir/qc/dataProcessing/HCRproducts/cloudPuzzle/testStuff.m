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