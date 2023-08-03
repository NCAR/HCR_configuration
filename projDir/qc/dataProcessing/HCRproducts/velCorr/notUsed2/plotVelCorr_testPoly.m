% Compare data from before and after velocity correction

clear all
close all

savefig=0;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc0'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

testCases=readtable(['/h/eol/romatsch/hcrCalib/velCorr/',project,'/testCases.dat']);

figdir=['/h/eol/romatsch/hcrCalib/velCorr/',project,'/velFigs/'];
formatOut = 'yyyymmdd_HHMM';

indir=HCRdir(project,quality,freqData);

polyTimePeriod=[2 5 10]; %Vector of time periods for poly fit in seconds
polyOrder=3; % Order of polynomial fit

plotWeatherLines=1; %Plot average velocity lines on velocity/alt plot

timeLag=0; %suspected time lag in seconds between HCR data and vertVel

for kk=1:size(testCases,1)
    startTime=datetime(testCases{kk,1},testCases{kk,2},testCases{kk,3}, ...
        testCases{kk,4},testCases{kk,5},0);
    endTime=datetime(testCases{kk,6},testCases{kk,7},testCases{kk,8}, ...
        testCases{kk,9},testCases{kk,10},0);
        
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.roll=[];
        data.pitch=[];
        data.tilt=[];
        data.rot=[];
        data.vertVel=[];
        data.elev=[];
        data.range=[];
        data.dbz=[];
        data.vel=[];
        data.velRaw=[];
        data.velCorr=[];
        
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
            data.roll=[data.roll,ncread(infile,'roll')'];
            data.pitch=[data.pitch,ncread(infile,'pitch')'];
            data.tilt=[data.tilt,ncread(infile,'tilt')'];
            data.rot=[data.rot,ncread(infile,'rotation')'];
            data.vertVel=[data.vertVel,ncread(infile,'vertical_velocity')'];
            
            data.dbz=[data.dbz,ncread(infile,'DBZ')];
            
            velIn=ncread(infile,'VEL');
            data.vel=[data.vel,velIn];
            velRawIn=ncread(infile,'VEL_RAW');
            data.velRaw=[data.velRaw,velRawIn];
            velCorrIn=ncread(infile,'VEL_CORR');
            data.velCorr=[data.velCorr,velCorrIn];
        end
        
        %% Get surface indices
        % Find surface indices
        DBZ=data.dbz;
        DBZ(1:15,:)=0;
        
        rangeTemp=data.range;
        DBZmask=DBZ; %Mask with DBZ values that are close to the surface
        
        rangeTemp(find(abs(data.range-data.alt)>150))=nan;
        DBZmask(isnan(rangeTemp))=nan;
        [bla ground_index]=nanmax(DBZmask,[],1);
        
        % Convert to linear indices
        linInd=sub2ind(size(DBZ),ground_index,1:length(data.alt));
        
        %% Adjust vel for time lag
        if timeLag>0
            vertVelUnlag=cat(2,nan(1,timeLag*10),data.vertVel(1:end-timeLag*10));
            data.vel=data.velRaw-vertVelUnlag;
        end
        %% Make poly fit
        velSmoothP=nan(length(polyTimePeriod),length(linInd));
        for ii=1:length(polyTimePeriod)
            velSmoothP(ii,:)=vel2vel_corr_testPoly(data.vel(linInd),data.time,polyTimePeriod(ii),polyOrder);
        end
        
        %% Average altitude bins
        %calculate above sea level altitudes
        asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
        
        maxAlt=max(max(data.alt))/1000;
        maxEdge=ceil( maxAlt/0.5 ) * 0.5;
        altEdges=0:0.5:maxEdge;
        
        velCorrP=nan(size(data.vel,1),size(data.vel,2),length(polyTimePeriod));
        for jj=1:length(polyTimePeriod)
            velCorrP(:,:,jj)=data.vel-velSmoothP(jj,:);
        end
        
        velRawMeans=nan((length(altEdges)-1),size(asl,2));
        velMeans=nan((length(altEdges)-1),size(asl,2));
        velCorrMeans=nan((length(altEdges)-1),size(asl,2));
        velCorrPMeans=nan((length(altEdges)-1),size(asl,2),length(polyTimePeriod));
        
        for ii=1:length(altEdges)-1
            binInd=find(asl>=altEdges(ii) & asl<altEdges(ii+1));
            
            velRawBin=nan(size(data.velRaw));
            velRawBin(binInd)=data.velRaw(binInd);
            velRawMeans(ii,:)=nanmean(velRawBin,1);
            
            velBin=nan(size(data.velRaw));
            velBin(binInd)=data.vel(binInd);
            velMeans(ii,:)=nanmean(velBin,1);
            
            velCorrBin=nan(size(data.velRaw));
            velCorrBin(binInd)=data.velCorr(binInd);
            velCorrMeans(ii,:)=nanmean(velCorrBin,1);
            
            for jj=1:length(polyTimePeriod)
                velCorrPBin=nan(size(data.velRaw));
                velCorrPpoly=velCorrP(:,:,jj);
                velCorrPBin(binInd)=velCorrPpoly(binInd);
                velCorrPMeans(ii,:,jj)=nanmean(velCorrPBin,1);
            end
        end
        
        middleLevel=round(maxEdge/2);
        %% Line plot
        
        close all
        f5=figure('DefaultAxesFontSize',12);
        set(f5,'Position',[200 500 3000 1500]);
        set(f5,'renderer','painters');
        
        subplot(3,1,1)
        hold on
        plot(data.time,data.vertVel,'-k','linewidth',1.5);
        plot(data.time,velRawMeans(middleLevel,:)-1,'-m','linewidth',1.5);
        plot(data.time,movmean(data.velRaw(linInd),50),'-g','linewidth',1.5);
        ylim([nanmedian(data.vertVel)-1.5 nanmedian(data.vertVel)+1.5]);
        ylabel('Velocity [m/s]');
        legend('VertVelPlane','VelRawWeather','VelRawSurfRunMean','Orientation','horizontal');
        title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut)],'interpreter','none');
        xlim([startTime,endTime]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        
        colorIn=lines(length(polyTimePeriod));
        legIn={'VelSurf','VelSurfRunMean'};
        
        subplot(3,1,2)
        hold on
        plot(data.time,data.vel(linInd),'color',[0.5 0.5 0.5],'linewidth',1.5);
        plot(data.time,movmean(data.vel(linInd),50),'-c','linewidth',1.5);
        for ii=1:length(polyTimePeriod)
            plot(data.time,velSmoothP(ii,:),'color',colorIn(ii,:),'linewidth',1.5);
            legIn=[legIn,num2str(polyTimePeriod(ii))];
        end
        plot(data.time,data.vel(linInd)-data.velCorr(linInd),'-m','linewidth',1.5);
        legIn=[legIn,'VelSurf-VelCorrSurf'];
        ylim([nanmedian(data.vel(linInd))-1 nanmedian(data.vel(linInd))+1]);
        ylabel('Velocity [m/s]');
        legend(legIn,'Orientation','horizontal');
        xlim([startTime,endTime]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        
        legIn2={'VelWeather'};
        
        subplot(3,1,3)
        hold on
        plot(data.time,velMeans(middleLevel,:),'-b','linewidth',1.5);
        for ii=1:length(polyTimePeriod)
            plot(data.time,velCorrPMeans(middleLevel,:,ii),'color',colorIn(ii,:),'linewidth',1.5);
            legIn2=[legIn2,num2str(polyTimePeriod(ii))];
        end
        legend(legIn2,'Orientation','horizontal');
        xlim([startTime,endTime]);
        ylim([nanmedian(velMeans(middleLevel,:))-1.5 nanmedian(velMeans(middleLevel,:))+1.5]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        ylabel('Velocity [m/s]');
        title('Weather');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f5, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_lines'],'-dpng','-r0');
        end
        
        %% Plot vel field
        
        if plotWeatherLines
            scaleFact=0.5;
            
            addVals=0.25:0.5:max(altEdges);
            
            velRawMadj=(velRawMeans-nanmean(velRawMeans,2))*scaleFact+addVals';
            velMadj=(velMeans-nanmean(velMeans,2))*scaleFact+addVals';
            velCorrMadj=(velCorrMeans-nanmean(velCorrMeans,2))*scaleFact+addVals';
        end
        
        f2=figure('DefaultAxesFontSize',12);
        set(f2,'Position',[200 500 1500 1500]);
        
        % Vel Raw
        subplot(3,1,1)
        hold on
        fig2=surf(data.time,asl,data.velRaw);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig2.CData;
        
        col_def1 = nan(size(fld));
        col_def2 = nan(size(fld));
        col_def3 = nan(size(fld));
        
        color_map=colormap(vel_default(29));
        %color_map=jet(17);
        
        %limits=[-8 -3.8 -3.3 -2.8 -2.3 -1.8 -1.3 -0.8 -0.3 0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 8];
        limits=-4.05:0.3:4.05;
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
        fig2.CData=col_def;
        
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
        
        if plotWeatherLines
            plot(data.time,velRawMadj,'-k');
        end
        title('VelRaw');
        
        % Vel
        
        subplot(3,1,2)
        hold on
        fig2=surf(data.time,asl,data.vel);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig2.CData;
        
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
        
        if plotWeatherLines
            plot(data.time,velMadj,'-k');
        end
        
        title('Vel');
        
        % VelCorr
        
        subplot(3,1,3)
        hold on
        fig2=surf(data.time,asl,data.velCorr);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig2.CData;
        
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
        
        if plotWeatherLines
            plot(data.time,velCorrMadj,'-k');
        end
        
        title('VelCorr');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_vel'],'-dpng','-r0');
        end
        
        %% Make poly fit plot
        if plotWeatherLines
            velCorrPMadj=(velCorrPMeans-nanmean(velCorrPMeans,2))*scaleFact+addVals';
        end
        
        f7=figure('DefaultAxesFontSize',12);
        set(f7,'Position',[1700 500 1500 1500]);
        
        for jj=1:length(polyTimePeriod)
            
            subplot(length(polyTimePeriod),1,jj)
            hold on
            fig4=surf(data.time,asl,velCorrP(:,:,jj));
            fig4.EdgeColor='none';
            ylim([-0.2 maxEdge]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig4.CData;
            
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
            fig4.CData=col_def;
            
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
            
            if plotWeatherLines
                plot(data.time,velCorrPMadj(:,:,jj),'-k');
            end
            
            title(['VelCorr poly fit 3, ',num2str(polyTimePeriod(jj)),' s']);
        end
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f7, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_velPoly'],'-dpng','-r0');
        end
    end
end