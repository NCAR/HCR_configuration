% Compare data from before and after velocity correction

clear all
close all

savefig=1;

flight='rf03';

startTime=datetime(2018,1,24,0,46,0);
endTime=datetime(2018,1,24,0,50,0);
%endTime=datetime(2018,1,16,0,21,0);

addpath('/h/eol/romatsch/git/private/utils/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
insdir='/scr/raf_data/SOCRATES/';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velFigs_weather/';
formatOut = 'yyyymmdd_HHMM';

fileList=makeFileList(indir,startTime,endTime,1);

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
        data.velDiff1=[data.velDiff1,velCorrIn-velIn];
    end
    
    velDiff11d=nanmax(data.velDiff1);
    velDiff11d(find(velDiff11d==0))=nan;
   
   
    %% Get surface data
    % Find surface indices
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
    
    velRawS=data.velRaw(linInd);
    velS=data.vel(linInd);
    
    %% Make poly fit
    polyTimePeriod=20; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP,velSmoothP]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    velCorrPS=velCorrP(linInd);
    
    %% Average altitude bins
    %calculate above sea level altitudes
    asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
    
    maxAlt=max(max(data.alt))/1000;
    maxEdge=ceil( maxAlt/0.5 ) * 0.5;
    altEdges=0:0.5:maxEdge;
    
    velRawMeans=nan((length(altEdges)-1),size(asl,2));
    velMeans=nan((length(altEdges)-1),size(asl,2));
    velCorrPMeans=nan((length(altEdges)-1),size(asl,2));
    velUnlagMeans=nan((length(altEdges)-1),size(asl,2));
    
    for ii=1:length(altEdges)-1
        binInd=find(asl>=altEdges(ii) & asl<altEdges(ii+1));
        
        velRawBin=nan(size(data.velRaw));
        velRawBin(binInd)=data.velRaw(binInd);
        velRawMeans(ii,:)=nanmean(velRawBin,1);
        
        velBin=nan(size(data.velRaw));
        velBin(binInd)=data.vel(binInd);
        velMeans(ii,:)=nanmean(velBin,1);
        
        velUnlagBin=nan(size(data.velRaw));
        velUnlagBin(binInd)=data.vel(binInd);
        velUnlagMeans(ii,:)=nanmean(velUnlagBin,1);
        
        velCorrPBin=nan(size(data.velRaw));
        velCorrPBin(binInd)=velCorrP(binInd);
        velCorrPMeans(ii,:)=nanmean(velCorrPBin,1);
    end
    
    %% Vert vel line and polyfit plot
    middleLevel=round(maxEdge/2);
    
    close all
    f5=figure('DefaultAxesFontSize',12);
    set(f5,'Position',[200 500 3000 500]);
    
    hold on
    plot(data.time,data.vertVel,'-r','linewidth',2);
    plot(data.time,velRawMeans(find(altEdges==middleLevel),:),'-b','linewidth',2);
    plot(data.time,velUnlagMeans(find(altEdges==middleLevel),:),'-k','linewidth',2);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(1):endTime];
    grid on
    grid minor
    ylabel('[m/s]');
    legend('VertVelPlatform','VelRaw 1-1.5km','Vel 1-1.5km','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' Mean vertical velocity between 1 and 1.5 km'],'interpreter','none');
    xlim([startTime,endTime]);
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f5, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_velRawToVel',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
    %%
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 3000 1000]);
    
    subplot(2,1,1)
    hold on
    plot(data.time,velS,'-c','linewidth',2);
    plot(data.time,-velCorrPS+velS,'-m','linewidth',2);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(1):endTime];
    grid on
    grid minor
    legend('VelOceanSurface','Filter','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' Vertical velocity of ocean surface'],'interpreter','none');
    xlim([startTime,endTime]);
    ylabel('[m/s]');
    
    subplot(2,1,2)
    hold on
    plot(data.time,velUnlagMeans(find(altEdges==middleLevel),:),'-k','linewidth',2);
    plot(data.time,velCorrPMeans(find(altEdges==middleLevel),:),'-g','linewidth',2);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(1):endTime];
    grid on
    grid minor
    legend('Vel 1-1.5km','VelCorr 1-1.5km','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' Mean vertical velocity between 3 and 3.5 km'],'interpreter','none');
    xlim([startTime,endTime]);
    ylabel('[m/s]');
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_velToVelCorr',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
     %%     
    close all
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 3000 500]);
    
    hold on
    plot(data.time,data.vertVel,'-r','linewidth',2);
    plot(data.time,velRawMeans(find(altEdges==middleLevel),:),'-b','linewidth',2);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(1):endTime];
    grid on
    grid minor
    ylabel('[m/s]');
    
    yyaxis right
    plot(data.time,data.alt,'-k','linewidth',2);
    ylabel('[m]');
    ylim([5642 5660])
    
    legend('VertVelPlatform','VelRaw 1-1.5km','AltitudePlatform','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' Mean vertical velocity between 1 and 1.5 km'],'interpreter','none');
    xlim([startTime,endTime]);
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_raf',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
end