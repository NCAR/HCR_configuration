% Remove backlobe echo

clear all
close all

savefig=1;

makefig4panels=0;
makefigMask=1;
makefigHist=0;
makefigCross=0;

whichCases='Weather';

%startTime=datetime(2018,2,20,4,17,0);
%endTime=datetime(2018,2,20,4,44,0);

addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/utils/colormaps/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
figdir=['/h/eol/romatsch/hcrCalib/backlobe/figs',whichCases,'/'];
formatOut = 'yyyymmdd_HHMM';

casesIn=readtable(['/h/eol/romatsch/hcrCalib/backlobe/infiles/cases',whichCases,'.txt']);
casesMat=table2array(casesIn);

% hist.dbz=[];
% hist.vel=[];
% hist.width=[];
% hist.ldr=[];
% edges.dbz=[];
% edges.width=[];
% edges.ldr=[];
% edges.vel=[];

for jj=1:size(casesIn,1)
    
    startTime=datetime(casesMat(jj,1:6));
    endTime=datetime(casesMat(jj,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.elev=[];
        data.range=[];
        data.dbz=[];
        data.velCorr=[];
        data.ldr=[];
        data.width=[];
        
        
        % Get uncorrected data
        for ii=1:size(fileList,2)
            infile=fileList{ii};
            
            startTimeIn=ncread(infile,'time_coverage_start')';
            startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
                str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
            timeRead=ncread(infile,'time')';
            timeIn=startTimeFile+seconds(timeRead);
            data.time=[data.time,timeIn];
            
            rangeIn=ncread(infile,'range');
            rangeMat=repmat(rangeIn,1,length(timeIn));
            data.range=[data.range,rangeMat];
            
            data.elev=[data.elev,ncread(infile,'elevation')'];
            data.alt=[data.alt,ncread(infile,'altitude')'];
            
            data.dbz=[data.dbz,ncread(infile,'DBZ')];
            data.ldr=[data.ldr,ncread(infile,'LDR')];
            data.velCorr=[data.velCorr,ncread(infile,'VEL_CORR')];
            data.width=[data.width,ncread(infile,'WIDTH_CORR')];
        end
        
        % Get right time frame
        startInd=find(data.time==startTime);
        endInd=find(data.time==endTime);
        
        data.time=data.time(:,startInd:endInd);
        data.alt=data.alt(:,startInd:endInd);
        data.elev=data.elev(:,startInd:endInd);
        data.range=data.range(:,startInd:endInd);
        data.dbz=data.dbz(:,startInd:endInd);
        data.velCorr=data.velCorr(:,startInd:endInd);
        data.ldr=data.ldr(:,startInd:endInd);
        data.width=data.width(:,startInd:endInd);
        
        %% Make backlobe echo mask
        
        % Initiate mask
        blMask=zeros(size(data.range)); % When using thresholds
        %blMask=ones(size(data.range)); % To test range
        
        % DBZ, VEL, LDR, SW threshold
        % blMask(data.dbz<-20 & data.velCorr<-4 & data.ldr>-12 & data.width>2)=1;
        blMask(data.dbz<-20 & data.width>1.5)=1;
        %blMask(data.width>3)=1;
        
        % Only within right altitude
        altMat=repmat(data.alt,size(data.range,1),1);
        % Lower limit
        blMask(data.range<(altMat-100))=0;
        % Upper limit
        blMask(data.range>(altMat+250))=0;
        
        % Only when scanning up
        blMask(:,data.elev<0)=0;
        
        %% Calculate above sea level altitudes
        asl=nan(size(data.range));
        downInd=find(data.elev<0);
        upInd=find(data.elev>=0);
        asl(:,downInd)=-1*((data.range(:,downInd).*cosd(abs(data.elev(downInd))-90)./1000)-data.alt(downInd)./1000);
        asl(:,upInd)=data.range(:,upInd).*cosd(abs(data.elev(upInd))-90)./1000+data.alt(upInd)./1000;
        
        %dbzInd=find(~isnan(data.dbz));
        %aslDBZ=asl(dbzInd);
        
        %maxAlt=max(max(data.alt))/1000;
        %maxEdge=ceil( maxAlt/0.5 ) * 0.5;
        
        %% Histogram data and plot
        close all
        maskInds=find(blMask==1);
        
        [dbzH dbzE ~]=histcounts(data.dbz(maskInds),50);
        [velH velE ~]=histcounts(data.velCorr(maskInds),50);
        [widthH widthE ~]=histcounts(data.width(maskInds),50);
        [ldrH ldrE ~]=histcounts(data.ldr(maskInds),50);
        
        dbzX=dbzE(1:end-1)+(dbzE(2)-dbzE(1))/2;
        velX=velE(1:end-1)+(velE(2)-velE(1))/2;
        widthX=widthE(1:end-1)+(widthE(2)-widthE(1))/2;
        ldrX=ldrE(1:end-1)+(ldrE(2)-ldrE(1))/2;
        
        if makefigHist
                        
            f3=figure('DefaultAxesFontSize',12);
            set(f3,'Position',[200 500 1000 1000]);
            
            % DBZ
            subplot(2,2,1)
            hold on
            bar(dbzX,dbzH,1)
            title('DBZ');
            
            % VEL
            subplot(2,2,2)
            hold on
            bar(velX,velH,1)
            title('VEL_CORR','interpreter','none');
            
            % WIDTH
            subplot(2,2,3)
            hold on
            bar(widthX,widthH,1)
            title('WIDTH_CORR','interpreter','none');
            
            % LDR
            subplot(2,2,4)
            hold on
            bar(ldrX,ldrH,1)
            title('LDR');
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f3, [figdir,'hist/',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_hist'],'-dpng','-r0');
            end
        end
        
        %% Cross section
        
        altRangeDiff=abs(altMat-data.range);
        
        [minDiff minRow]=min(altRangeDiff,[],1);
        
        rowPlus=15;
        rowMinus=6;
        
        crossData.dbz=nan(rowPlus+rowMinus+1,length(minRow));
        crossData.vel=nan(rowPlus+rowMinus+1,length(minRow));
        crossData.width=nan(rowPlus+rowMinus+1,length(minRow));
        crossData.ldr=nan(rowPlus+rowMinus+1,length(minRow));
        
        for kk=1:length(minRow)
            crossData.dbz(:,kk)=data.dbz(minRow(kk)-rowMinus:minRow(kk)+rowPlus,kk);
            crossData.vel(:,kk)=data.velCorr(minRow(kk)-rowMinus:minRow(kk)+rowPlus,kk);
            crossData.width(:,kk)=data.width(minRow(kk)-rowMinus:minRow(kk)+rowPlus,kk);
            crossData.ldr(:,kk)=data.ldr(minRow(kk)-rowMinus:minRow(kk)+rowPlus,kk);
        end
        
        crossData.dbz(:,data.elev<0)=nan;
        crossData.vel(:,data.elev<0)=nan;
        crossData.width(:,data.elev<0)=nan;
        crossData.ldr(:,data.elev<0)=nan;
        
        ind.dbz=zeros(size(crossData.dbz));
        ind.dbz(~isnan(crossData.dbz))=1;
        ind.vel=zeros(size(crossData.dbz));
        ind.vel(~isnan(crossData.vel))=1;
        ind.width=zeros(size(crossData.dbz));
        ind.width(~isnan(crossData.width))=1;
        ind.ldr=zeros(size(crossData.dbz));
        ind.ldr(~isnan(crossData.ldr))=1;
        
        if makefigCross
                        
            f4=figure('DefaultAxesFontSize',12);
            set(f4,'Position',[200 500 1500 1000]);
            
            % DBZ
            subplot(2,4,1)
            hold on
            bar(nanmean(crossData.dbz,2),1)
            title('DBZ mean');
            xlim([1 size(ind.dbz,1)])
            
            subplot(2,4,2)
            hold on
            bar(nansum(ind.dbz,2),1)
            title('DBZ pixels');
            xlim([1 size(ind.dbz,1)])
            
            % VEL
            subplot(2,4,3)
            hold on
            bar(nanmean(crossData.vel,2),1)
            title('VEL_CORR mean','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            subplot(2,4,4)
            hold on
            bar(nansum(ind.vel,2),1)
            title('VEL_CORR pixels','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            % WIDTH
            subplot(2,4,5)
            hold on
            bar(nanmean(crossData.width,2),1)
            title('WIDTH_CORR mean','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            subplot(2,4,6)
            hold on
            bar(nansum(ind.width,2),1)
            title('WIDTH_CORR pixels','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            % LDR
            subplot(2,4,7)
            hold on
            bar(nanmean(crossData.ldr,2),1)
            title('LDR mean','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            subplot(2,4,8)
            hold on
            bar(nansum(ind.ldr,2),1)
            title('LDR pixels','interpreter','none');
            xlim([1 size(ind.dbz,1)])
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f4, [figdir,'cross/',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_cross'],'-dpng','-r0');
            end
        end
        %% Plot 4 panels with variables
        
        if makefig4panels
                        
            f1=figure('DefaultAxesFontSize',12);
            set(f1,'Position',[200 500 2000 1500]);
            
            % DBZ
            subplot(2,2,1)
            hold on
            
            fig1=surf(data.time,asl,data.dbz);
            fig1.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig1.CData;
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            limits=[-inf (-43:3:23) inf];
            color_map=colormap(dbz_default);
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
            fig1.CData=col_def;
            
            hcb=colorbar;
            set(get(hcb,'Title'),'String','dBZ');
            colormap(gca,color_map);
            caxis([0 size(color_map,1)]);
            caxis_yticks=(1:1:size(color_map,1)-1);
            caxis_ytick_labels=num2str(limits(2:end-1)');
            while length(caxis_yticks)>16
                caxis_yticks=caxis_yticks(1:2:end);
                caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
            end
            set(hcb,'ytick',caxis_yticks);
            set(hcb,'YTickLabel',caxis_ytick_labels);
            ylabel('Altitude [km]');
            
            title('DBZ');
            
            % VEL_CORR
            subplot(2,2,2)
            hold on
            
            fig1=surf(data.time,asl,data.velCorr);
            fig1.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig1.CData;
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            color_map=colormap(vel_default(18));
            limits=[-8 -3.8 -3.3 -2.8 -2.3 -1.8 -1.3 -0.8 -0.3 0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 8];
            limits=[-inf limits inf];
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
            fig1.CData=col_def;
            
            hcb=colorbar;
            set(get(hcb,'Title'),'String','m/s');
            colormap(gca,color_map);
            caxis([0 size(color_map,1)]);
            caxis_yticks=(1:1:size(color_map,1)-1);
            caxis_ytick_labels=num2str(limits(2:end-1)');
            while length(caxis_yticks)>16
                caxis_yticks=caxis_yticks(1:2:end);
                caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
            end
            set(hcb,'ytick',caxis_yticks);
            set(hcb,'YTickLabel',caxis_ytick_labels);
            ylabel('Altitude [km]');
            
            title('VEL_CORR','interpreter','none');
            
            % WIDTH
            subplot(2,2,3)
            hold on
            
            fig1=surf(data.time,asl,data.width);
            fig1.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig1.CData;
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            color_map=colormap(width_default);
            limits=[0.1 0.2 0.3 0.4 0.5 0.6 0.75 1 1.25 1.5 1.75 2 2.25 2.5 3 4];
            limits=[-inf limits inf];
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
            fig1.CData=col_def;
            
            hcb=colorbar;
            set(get(hcb,'Title'),'String','m/s');
            colormap(gca,color_map);
            caxis([0 size(color_map,1)]);
            caxis_yticks=(1:1:size(color_map,1)-1);
            caxis_ytick_labels=num2str(limits(2:end-1)');
            while length(caxis_yticks)>16
                caxis_yticks=caxis_yticks(1:2:end);
                caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
            end
            set(hcb,'ytick',caxis_yticks);
            set(hcb,'YTickLabel',caxis_ytick_labels);
            ylabel('Altitude [km]');
            
            title('WIDTH_CORR','interpreter','none');
            
            % LDR
            subplot(2,2,4)
            hold on
            
            fig1=surf(data.time,asl,data.ldr);
            fig1.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig1.CData;
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            %color_map=parula(25);
            color_map=colormap(width_default(25));
            %color_map=cat(1,[0.2 0.2 0.2],color_map);
            color_map=flipud(color_map);
            limits=[-35:1:-23,-22:2:-10,-6:4:6];
            limits=[-inf limits inf];
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
            fig1.CData=col_def;
            
            hcb=colorbar;
            set(get(hcb,'Title'),'String','dB');
            colormap(gca,color_map);
            caxis([0 size(color_map,1)]);
            caxis_yticks=(1:1:size(color_map,1)-1);
            caxis_ytick_labels=num2str(limits(2:end-1)');
            while length(caxis_yticks)>16
                caxis_yticks=caxis_yticks(1:2:end);
                caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
            end
            set(hcb,'ytick',caxis_yticks);
            set(hcb,'YTickLabel',caxis_ytick_labels);
            ylabel('Altitude [km]');
            
            title('LDR');
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f1, [figdir,'vars/',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_vars'],'-dpng','-r0');
            end
        end
        
        %% Mask figure
        
        if makefigMask
            
            f2=figure('DefaultAxesFontSize',12);
            set(f2,'Position',[200 500 1500 1500]);
            
            subplot(3,1,1)
            hold on
            fig2=surf(data.time,asl,data.dbz);
            fig2.EdgeColor='none';
            fig2.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig2.CData;
            
            limits=[-inf (-43:3:23) inf];
            color_map=colormap(dbz_default);
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
                        
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
                        
            fig2.CData=col_def;
            
           hcb=colorbar;
           set(get(hcb,'Title'),'String','dBZ');
           colormap(gca,color_map);
           caxis([0 size(color_map,1)]);
           caxis_yticks=(1:1:size(color_map,1)-1);
           caxis_ytick_labels=num2str(limits(2:end-1)');
          while length(caxis_yticks)>16
               caxis_yticks=caxis_yticks(1:2:end);
               caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
           end
           set(hcb,'ytick',caxis_yticks);
           set(hcb,'YTickLabel',caxis_ytick_labels);
            
            ylabel('Altitude [km]');
            
            title('All data');
            
            subplot(3,1,2)
            hold on
            fig2=surf(data.time,asl,data.dbz);
            fig2.EdgeColor='none';
            fig2.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig2.CData;
            
            limits=[-inf (-43:3:23) inf];
            color_map=colormap(dbz_default);
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            col_def1(blMask==1)=1;
            col_def2(blMask==1)=0;
            col_def3(blMask==1)=0;
            
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
                        
            fig2.CData=col_def;
            
           hcb=colorbar;
           set(get(hcb,'Title'),'String','dBZ');
           colormap(gca,color_map);
           caxis([0 size(color_map,1)]);
           caxis_yticks=(1:1:size(color_map,1)-1);
           caxis_ytick_labels=num2str(limits(2:end-1)');
          while length(caxis_yticks)>16
               caxis_yticks=caxis_yticks(1:2:end);
               caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
           end
           set(hcb,'ytick',caxis_yticks);
           set(hcb,'YTickLabel',caxis_ytick_labels);
            
            ylabel('Altitude [km]');
            
            title('With mask');
            
            subplot(3,1,3)
            hold on
            fig2=surf(data.time,asl,data.dbz);
            fig2.EdgeColor='none';
            fig2.EdgeColor='none';
            ylim([-1 3.5]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig2.CData;
            
            limits=[-inf (-43:3:23) inf];
            color_map=colormap(dbz_default);
            
            col_def1 = nan(size(fld));
            col_def2 = nan(size(fld));
            col_def3 = nan(size(fld));
            
            for ii=1:size(color_map,1)
                col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
                col_def1(col_ind)=color_map(ii,1);
                col_def2(col_ind)=color_map(ii,2);
                col_def3(col_ind)=color_map(ii,3);
            end
            col_def1(blMask==1)=1;
            col_def2(blMask==1)=1;
            col_def3(blMask==1)=1;
            
            if ~isequal(size(col_def1),(size(fld)))
                col_def=cat(3,col_def1',col_def2',col_def3');
            else
                col_def=cat(3,col_def1,col_def2,col_def3);
            end
                        
            fig2.CData=col_def;
            
           hcb=colorbar;
           set(get(hcb,'Title'),'String','dBZ');
           colormap(gca,color_map);
           caxis([0 size(color_map,1)]);
           caxis_yticks=(1:1:size(color_map,1)-1);
           caxis_ytick_labels=num2str(limits(2:end-1)');
          while length(caxis_yticks)>16
               caxis_yticks=caxis_yticks(1:2:end);
               caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
           end
           set(hcb,'ytick',caxis_yticks);
           set(hcb,'YTickLabel',caxis_ytick_labels);
            
            ylabel('Altitude [km]');
            
            title('Masked out');
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f2, [figdir,'maskWidthThresh/',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_mask'],'-dpng','-r0');
            end
        end
    end
end