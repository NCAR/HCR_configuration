% Compare data from before and after velocity correction

clear all
close all

savefig=0;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

%testCases=readtable('/h/eol/romatsch/hcrCalib/velCorr/testCases.dat');
testCases=readtable(['/h/eol/romatsch/hcrCalib/velCorr/otrec/testCases.dat']);

for kk=1:size(testCases,1)
    startTime=datetime(testCases{kk,1},testCases{kk,2},testCases{kk,3}, ...
        testCases{kk,4},testCases{kk,5},0);
    endTime=datetime(testCases{kk,6},testCases{kk,7},testCases{kk,8}, ...
        testCases{kk,9},testCases{kk,10},0);
    
    % startTime=datetime(2018,1,16,0,20,0);
    % endTime=datetime(2018,1,16,0,25,0);
        
    indir='/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/';
    %figdir='/h/eol/romatsch/hcrCalib/velCorr/velDiffFigs_compare_cases/';
    %figdir='/h/eol/romatsch/hcrCalib/velCorr/velFigs_wholeProcess2/';
    figdir=['/h/eol/romatsch/hcrCalib/velCorr/otrec/velFigs/scott/'];
    formatOut = 'yyyymmdd_HHMM';
    
%     surfVelFile='/scr/rain1/rsfdata/projects/socrates/hcr/qc/txt/HcrVelCorrect.txt';
%     surfVelAll=txtTable2matTable(surfVelFile,' ');
%     surfVelTime=datetime(surfVelAll.year,surfVelAll.month,surfVelAll.day,...
%         surfVelAll.hour,surfVelAll.min,surfVelAll.sec);
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.elev=[];
        data.range=[];
        data.dbz=[];
        %data.snrvc=[];
        data.vel=[];
        data.velCorr=[];
        data.velDiff1=[];
        
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
            %data.snrvc=[data.snrvc,ncread(infile,'SNRVC')];
            
            velIn=ncread(infile,'VEL');
            data.vel=[data.vel,velIn];
            %velCorrIn=ncread(infile,'VEL_CORR');
%            data.velCorr=[data.velCorr,velCorrIn];
  %          data.velDiff1=[data.velDiff1,velCorrIn-velIn];
        end
        
        velDiff11d=nanmax(data.velDiff1);
        velDiff11d(find(velDiff11d==0))=nan;
        
        % Surface velocity data
%        surfVelInds=find(surfVelTime>=startTime & surfVelTime<=endTime);
        
        %% Get Scott's data
        [velCorrS,velCorrSurfS,rangeToSurfS,velSmoothS]=vel2vel_corr(data.vel,data.dbz,data.range,data.alt);
        
        % Plot line plots
        
        close all
        
        %     f1=figure('DefaultAxesFontSize',14);
        %     set(f1,'Position',[200 500 1500 1500]);
        %
        %     subplot(4,1,1)
        %     hold on
        %     plot(data.time,velDiff11d,'-r','linewidth',1.2);
        %     plot(data.time,movmean(velDiff11d,500),'-k','linewidth',3);
        %     ylim([-2 2]);
        %     ylabel('VelCorr - Vel [m/s]');
        %     legend('VelCorr - Vel (Mike)','Orientation','horizontal');
        %     title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut)],'interpreter','none');
        %     xlim([startTime,endTime]);
        %     grid on
        %
        %     if ~isempty(surfVelInds)
        %         subplot(4,1,2)
        %         hold on
        %         plot(surfVelTime(surfVelInds),surfVelAll.VelMeas(surfVelInds),'-b');
        %         plot(surfVelTime(surfVelInds),surfVelAll.VelFilt(surfVelInds),'-r','linewidth',2);
        %         plot(data.time,velSmoothS,'-k','linewidth',2);
        %         ylim([-2 2]);
        %         xlim([startTime,endTime]);
        %         ylabel('Surface velocity [m/s]');
        %         yticks([-2 -1 0 1 2]);
        %         legend('Vel','VelSmooth Mike','VelSmooth Scott','Orientation','horizontal');
        %         grid on
        %
        %         DBZ=data.dbz;
        %         DBZ(1:15,:)=0;
        %
        %         rangeTemp=data.range;
        %         DBZmask=DBZ; %Mask with DBZ values that are close to the surface
        %
        %         rangeTemp(find(abs(data.range-data.alt)>150))=nan;
        %         DBZmask(isnan(rangeTemp))=nan;
        %         [bla ground_index]=nanmax(DBZmask,[],1);
        %         wrong_ground_ind=find(ground_index==1);
        %
        %         % Convert to linear indices
        %         linInd=sub2ind(size(DBZ),ground_index,1:length(data.alt));
        %
        %         surfVelData=data.velCorr(linInd);
        %         surfVelData(wrong_ground_ind)=nan;
        %
        %         subplot(4,1,3)
        %         hold on
        %         plot(data.time,velCorrSurfS,'-r','linewidth',1);
        %         plot(data.time,surfVelData,'-b','linewidth',1);
        %         plot(surfVelTime(surfVelInds),surfVelAll.VelCorr(surfVelInds),'-k','linewidth',1);
        %         %plot(data.time,velCorrS,'-r','linewidth',1);
        %         ylim([-2 2]);
        %         ylabel('Surface velocity [m/s]');
        %         yticks([-2 -1 0 1 2]);
        %         xlim([startTime,endTime]);
        %         %legend('VelSurf Mike','VelSurf Scott','Orientation','horizontal');
        %         legend('VelSurf Scott','VelSurf Data','VelSurf Mike','Orientation','horizontal');
        %         grid on
        %     end
        %
        %     subplot(4,1,4);
        %     hold on;
        %     plot(data.time,data.elev);
        %     ylabel('Elev angle [deg]');
        %     ylim([-100 210]);
        %     xlim([startTime,endTime]);
        %     yticks([-90 -45 0 45 90 135 180]);
        %     grid on
        %
        %     yyaxis right
        %     plot(data.time,data.alt./1000);
        %     if ~isempty(surfVelInds)
        %         plot(surfVelTime(surfVelInds),surfVelAll.RangeToSurface(surfVelInds),'-c');
        %     end
        %     plot(data.time,rangeToSurfS./1000,'-k');
        %     ylim([0 9]);
        %     ylabel('Altitude/Range [km]');
        %     yticks([0 1.5 3 4.5 6 7.5 9]);
        %     legend('ElevAng','Alt','RangeToSurf Mike','RangeToSurf Scott','Orientation','horizontal');
        %
        %     if savefig
        %         set(gcf,'PaperPositionMode','auto')
        %         print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
        %     end
        
        % Plot vel field
        
        %calculate above sea level altitudes
        asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
        
        maxAlt=max(max(data.alt))/1000;
        maxEdge=ceil( maxAlt/0.5 ) * 0.5;
        
        f2=figure('DefaultAxesFontSize',14);
        set(f2,'Position',[200 500 1500 1500]);
        
        subplot(3,1,1)
        fig2=surf(data.time,asl,data.vel);
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
        title('Vel');
        
        subplot(3,1,2)
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
        title('VelCorr Mike');
        
        
        subplot(3,1,3)
        fig3=surf(data.time,asl,velCorrS);
        fig3.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig3.CData;
        
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
        fig3.CData=col_def;
        
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
        title('VelCorr Scott');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_velCompare'],'-dpng','-r0');
        end
        
    end
end