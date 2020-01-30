% Compare data from before and after velocity correction

clear all
close all

startTimeFlight=datetime(2018,1,16,0,0,0); %Whole hour in which flight starts

%flightHours=1;
flightHours=24;

addpath('/h/eol/romatsch/git/private/utils/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velDiffFigs_compare/';
formatOut = 'yyyymmdd_HH';

surfVelFile='/scr/rain1/rsfdata/projects/socrates/hcr/qc/txt/HcrVelCorrect.txt';
surfVelAll=txtTable2matTable(surfVelFile,' ');
surfVelTime=datetime(surfVelAll.year,surfVelAll.month,surfVelAll.day,...
    surfVelAll.hour,surfVelAll.min,surfVelAll.sec);

for jj=1:flightHours
    
    startTime=startTimeFlight+hours(jj-1);
    endTime=startTimeFlight+hours(jj);
    
    fileList=makeFileList(indir,startTime,endTime,1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.elev=[];
        data.range=[];
        data.dbz=[];
        data.snrvc=[];
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
            data.snrvc=[data.snrvc,ncread(infile,'SNRVC')];
            
            velIn=ncread(infile,'VEL');
            data.vel=[data.vel,velIn];
            velCorrIn=ncread(infile,'VEL_CORR');
            data.velCorr=[data.velCorr,velCorrIn];
            data.velDiff1=[data.velDiff1,velCorrIn-velIn];
        end
        
        velDiff11d=nanmax(data.velDiff1);
        velDiff11d(find(velDiff11d==0))=nan;
       
        % Surface velocity data
        surfVelInds=find(surfVelTime>=startTime & surfVelTime<=endTime);
        
        %% Get Scott's data
        [velCorrS,velCorrSurfS,rangeToSurfS]=vel2vel_corr(data.vel,data.dbz,data.range,data.alt);
        
        %% Plot        
                
        close all
        
        f1=figure('DefaultAxesFontSize',14);
        set(f1,'Position',[200 500 1500 1500]);
        
        subplot(4,1,1)
        hold on
        plot(data.time,velDiff11d,'-r','linewidth',1.2);
        plot(data.time,movmean(velDiff11d,500),'-k','linewidth',3);
        ylim([-2 2]);
        ylabel('VelCorr - Vel [m/s]');
        legend('VelCorr - Vel (Mike)','Orientation','horizontal');
        title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut)],'interpreter','none');
        xlim([startTime,endTime]);
        grid on
        
        if ~isempty(surfVelInds)
            subplot(4,1,2)
            hold on
            plot(surfVelTime(surfVelInds),surfVelAll.VelMeas(surfVelInds),'-b');
            plot(surfVelTime(surfVelInds),surfVelAll.VelFilt(surfVelInds),'-r');
            plot(surfVelTime(surfVelInds),surfVelAll.VelCorr(surfVelInds),'-k','linewidth',2);
            ylim([-2 2]);
            xlim([startTime,endTime]);
            ylabel('Surface velocity [m/s]');
            yticks([-2 -1 0 1 2]);
            legend('Vel','VelSmooth','VelCorr','Orientation','horizontal');
            grid on
            
            DBZ=data.dbz;
            DBZ(1:15,:)=0;
            
            rangeTemp=data.range;
            DBZmask=DBZ; %Mask with DBZ values that are close to the surface
            
            rangeTemp(find(abs(data.range-data.alt)>150))=nan;
            DBZmask(isnan(rangeTemp))=nan;
            [bla ground_index]=nanmax(DBZmask,[],1);
            wrong_ground_ind=find(ground_index==1);
            
            % Convert to linear indices
            linInd=sub2ind(size(DBZ),ground_index,1:length(data.alt));
            
            surfVelData=data.velCorr(linInd);
            surfVelData(wrong_ground_ind)=nan;
            
            subplot(4,1,3)
            hold on
            plot(data.time,velCorrSurfS,'-r','linewidth',1);
            plot(data.time,surfVelData,'-b','linewidth',1);
            plot(surfVelTime(surfVelInds),surfVelAll.VelCorr(surfVelInds),'-k','linewidth',1);
            %plot(data.time,velCorrS,'-r','linewidth',1);
            ylim([-2 2]);
            ylabel('Surface velocity [m/s]');
            yticks([-2 -1 0 1 2]);
            xlim([startTime,endTime]);
            legend('VelSurf Scott','VelSurf data','VelSurf Mike','Orientation','horizontal');
            grid on
        end
        
        subplot(4,1,4);
        hold on;
        plot(data.time,data.elev);
        ylabel('Elev angle [deg]');
        ylim([-100 210]);
        xlim([startTime,endTime]);
        yticks([-90 -45 0 45 90 135 180]);
        grid on
        
        yyaxis right
        plot(data.time,data.alt./1000);
        if ~isempty(surfVelInds)
            plot(surfVelTime(surfVelInds),surfVelAll.RangeToSurface(surfVelInds),'-c');
        end
        plot(data.time,rangeToSurfS./1000,'-k');
        ylim([-0.23 6.68]);
        ylabel('Altitude/Range [km]');
        yticks([0 1 2 3 4 5 6]);
        legend('ElevAng','Alt','RangeToSurf Mike','RangeToSurf Scott','Orientation','horizontal');
        
        set(gcf,'PaperPositionMode','auto')
        print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
        
    end
end