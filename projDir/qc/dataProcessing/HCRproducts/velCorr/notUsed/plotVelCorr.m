% Compare data from before and after velocity correction

clear all
close all

startTimeFlight=datetime(2018,1,15,0,0,0); %Whole hour in which flight starts

flightHours=1000;

addpath('/h/eol/romatsch/git/private/utils/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velDiffFigs/';
formatOut = 'yyyymmdd_HH';

for jj=1:flightHours
    
    startTime=startTimeFlight+hours(jj-1);
    endTime=startTimeFlight+hours(jj);
    
    fileList=makeFileList(indir,startTime,endTime,1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.range=[];
        data.elev=[];
        %data.az=[];
        data.roll=[];
        data.pitch=[];
        data.drift=[];
        data.rot=[];
        data.tilt=[];
        data.velRaw=[];
        data.vel=[];
        data.velCorr=[];
        data.velDiff1=[];
        data.velDiff2=[];
        data.velDiff3=[];
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
            %data.az=[data.az,ncread(infile,'azimuth')'];
            data.roll=[data.roll,ncread(infile,'roll')'];
            data.pitch=[data.pitch,ncread(infile,'pitch')'];
            data.drift=[data.drift,ncread(infile,'drift')'];
            data.rot=[data.rot,ncread(infile,'rotation')'];
            data.tilt=[data.tilt,ncread(infile,'tilt')'];
            
            data.alt=[data.alt,ncread(infile,'altitude')'];
            
            velIn=ncread(infile,'VEL');
            data.vel=[data.vel,velIn];
            velRawIn=ncread(infile,'VEL_RAW');
            data.velRaw=[data.velRaw,velRawIn];
            velCorrIn=ncread(infile,'VEL_CORR');
            data.velCorr=[data.velCorr,velCorrIn];
            data.velDiff1=[data.velDiff1,velCorrIn-velIn];
            data.velDiff2=[data.velDiff2,velCorrIn-velRawIn];
            data.velDiff3=[data.velDiff3,velIn-velRawIn];
        end
        
        velDiff11d=nanmax(data.velDiff1);
        velDiff11d(find(velDiff11d==0))=nan;
        
        velDiff21d=nanmax(data.velDiff2);
        velDiff21d(find(velDiff21d==0))=nan;
        
        velDiff31d=nanmax(data.velDiff3);
        velDiff31d(find(velDiff31d==0))=nan;
                        
        close all
        
        f1=figure('DefaultAxesFontSize',14);
        set(f1,'Position',[200 500 1500 1500]);
        
        subplot(3,1,1)
        hold on
        plot(data.time,velDiff21d,'-b');
        plot(data.time,velDiff31d,'-c');
        ylim([-12 12]);
        ylabel('Velocity correction [m/s]');
        yticks([-12 -9 -6 -3 0 3 6 9 12]);
        
        yyaxis right
        plot(data.time,velDiff11d,'-r','linewidth',1.2);
        plot(data.time,movmean(velDiff11d,500),'-k','linewidth',3);
        ylim([-2 2]);
        ylabel('VelCorr - Vel [m/s]');
        ax = gca;
        ax.YColor=[1 0 0];
        legend('VelCorr - VelRaw','Vel - VelRaw','VelCorr - Vel','Orientation','horizontal');
        title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut)],'interpreter','none');
        xlim([startTime,endTime]);
        grid on
        
        subplot(3,1,2);
        hold on;
        plot(data.time,data.elev);
        rotPlot=data.rot;
        rotPlot(find(rotPlot>350))=rotPlot(find(rotPlot>350))-360;
        plot(data.time,rotPlot);
        ylabel('Elev and rot [deg]');
        ylim([-100 210]);
        xlim([startTime,endTime]);
        yticks([-90 -45 0 45 90 135 180]);
        grid on
        
        yyaxis right
        plot(data.time,abs(data.elev)-90);
        ylabel('Abs(Elev)-90 [deg]');
        ylim([-0.485 0.545])
        yticks([-0.45 -0.3 -0.15 0 0.15 0.3 0.45]);
        
        legend('Elev','Rot','Abs(Elev)-90','Orientation','horizontal');
        
        subplot(3,1,3);
        hold on;
        plot(data.time,data.pitch);
        plot(data.time,data.tilt);
        plot(data.time,data.roll);
        plot(data.time,data.drift,'-g');
        ylabel('Angles [deg]');
        ylim([-30 30]);
        xlim([startTime,endTime]);
        grid on
        
        yyaxis right
        plot(data.time,data.alt./1000);
        ylim([0 9]);
        ylabel('Altitude [km]');
        yticks([0 1.5 3 4.5 6 7.5 9]);
        legend('Pitch','Tilt','Roll','Drift','Alt','Orientation','horizontal');
        
        set(gcf,'PaperPositionMode','auto')
        print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
    end
end