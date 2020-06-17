% Compare data from before and after velocity correction

clear all
close all

savefig=0;

flight='rf03';

startTime=datetime(2018,1,24,0,45,0);
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
        
    %% Get INS data
    insfile=[insdir,'SOCRATES',flight,'srt.nc'];
    info=ncinfo(insfile);
    
    insstartTimeIn=info.Variables(1).Attributes(3).Value;
    insstartTimeFile=datetime(str2num(insstartTimeIn(15:18)),str2num(insstartTimeIn(20:21)),str2num(insstartTimeIn(23:24)),...
        str2num(insstartTimeIn(26:27)),str2num(insstartTimeIn(29:30)),str2num(insstartTimeIn(32:33)));
    instimeRead=ncread(insfile,'Time')';
    instime=insstartTimeFile+seconds(instimeRead);
    
    insInds=find(instime>=startTime & instime<=endTime);
    instime=instime(insInds);
    
    % 50 hz time
    sps50add=0:1/50:1;
    sps50add=sps50add(1:end-1);
    
    timeMat50=repmat(instime,50,1);
    timeMat50add=timeMat50+repmat(seconds(sps50add'),1,size(timeMat50,2));
    
    instimeVec50=reshape(timeMat50add,1,[]);
    
    % 10 hz time
    sps10add=0:1/10:1;
    sps10add=sps10add(1:end-1);
    
    timeMat10=repmat(instime,10,1);
    timeMat10add=timeMat10+repmat(seconds(sps10add'),1,size(timeMat10,2));
    
    instimeVec10=reshape(timeMat10add,1,[]);
       
    irs1=ncread(insfile,'VSPD');
    irs1=irs1(:,insInds);
    irs2=ncread(insfile,'VSPD_IRS2');
    irs2=irs2(:,insInds);
    irs3=ncread(insfile,'VSPD_IRS3');
    irs3=irs3(:,insInds);
    gps1=ncread(insfile,'GGVSPD');
    gps1=gps1(:,insInds);
    gps2=ncread(insfile,'GGVSPD_GPS1');
    gps2=gps2(:,insInds);
    
    irs1Vec=reshape(irs1,1,[]);
    irs2Vec=reshape(irs2,1,[]);
    irs3Vec=reshape(irs3,1,[]);
    gps1Vec=reshape(gps1,1,[]);
    gps2Vec=reshape(gps2,1,[]);
   
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
    
    %% Unlag vertVel
    timeLag=0; %suspected time lag in seconds between HCR data and vertVel
    vertVelUnlag=cat(2,nan(1,timeLag*10),data.vertVel(1:end-timeLag*10));
    data.velUnlag=data.velRaw-vertVelUnlag;
    velUnlagS=data.velUnlag(linInd);
    %% Make poly fit
    polyTimePeriod=20; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP,velSmoothP]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    velCorrPS=velCorrP(linInd);
    
    polyTimePeriod=30; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP30,velSmoothP30]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    velCorrP30S=velCorrP30(linInd);
    
    polyTimePeriod=60; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP60,velSmoothP60]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    velCorrP60S=velCorrP60(linInd);
    
    polyTimePeriod=120; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP120,velSmoothP120]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    velCorrP120S=velCorrP120(linInd);
    %% Average altitude bins
    %calculate above sea level altitudes
    asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
    
    maxAlt=max(max(data.alt))/1000;
    maxEdge=ceil( maxAlt/0.5 ) * 0.5;
    altEdges=0:0.5:maxEdge;
    
    velRawMeans=nan((length(altEdges)-1),size(asl,2));
    velMeans=nan((length(altEdges)-1),size(asl,2));
    velCorrPMeans=nan((length(altEdges)-1),size(asl,2));
    velCorrP30Means=nan((length(altEdges)-1),size(asl,2));
    velCorrP60Means=nan((length(altEdges)-1),size(asl,2));
    velCorrP120Means=nan((length(altEdges)-1),size(asl,2));
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
        velUnlagBin(binInd)=data.velUnlag(binInd);
        velUnlagMeans(ii,:)=nanmean(velUnlagBin,1);
        
        velCorrPBin=nan(size(data.velRaw));
        velCorrPBin(binInd)=velCorrP(binInd);
        velCorrPMeans(ii,:)=nanmean(velCorrPBin,1);
        
        velCorrP30Bin=nan(size(data.velRaw));
        velCorrP30Bin(binInd)=velCorrP30(binInd);
        velCorrP30Means(ii,:)=nanmean(velCorrP30Bin,1);
        
        velCorrP60Bin=nan(size(data.velRaw));
        velCorrP60Bin(binInd)=velCorrP60(binInd);
        velCorrP60Means(ii,:)=nanmean(velCorrP60Bin,1);
        
        velCorrP120Bin=nan(size(data.velRaw));
        velCorrP120Bin(binInd)=velCorrP120(binInd);
        velCorrP120Means(ii,:)=nanmean(velCorrP120Bin,1);
    end
    
    %% Vert vel line and polyfit plot
    middleLevel=round(maxEdge/2);
    
    %close all
    f5=figure('DefaultAxesFontSize',12);
    set(f5,'Position',[200 500 3000 1500]);
    
    subplot(3,1,1)
    hold on
    plot(instimeVec50,irs1Vec,'-b');
    plot(instimeVec50,irs2Vec,'-c');
    plot(instimeVec50,irs3Vec,'-g');
    plot(instimeVec10,gps1Vec,'-k');
    plot(instimeVec10,gps2Vec,'-','color',[0.5 0.5 0.5]);
    plot(data.time,data.vertVel,'-r');
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(1):endTime];
    grid on
    grid minor
    legend('IRS1','IRS2','IRS3','GPS1','GPS2','VertVel','Orientation','horizontal');
    
    subplot(3,1,2)
    hold on
    plot(data.time,velUnlagS,'color',[0.5 0.5 0.5]);
    plot(data.time,movmean(velUnlagS,50),'color',[1 0 0]);
    plot(data.time,velSmoothP,'-k');
    ylim([-2 2]);
    ylabel('Vel smooth [m/s]');
    legend('VelSurf','VelSurfRunMean','VelSmooth Poly3','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' poly fit ',num2str(polyTimePeriod),' s'],'interpreter','none');
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
    grid on
    grid minor
    
    subplot(3,1,3)
    hold on
    plot(data.time,velRawMeans(middleLevel,:)-1,'-m');
    plot(data.time,data.vertVel,'-c','linewidth',2);
    plot(data.time,movmean(velRawS,50),'-r');
    ylim([-2 2]);
    ylabel('Vel smooth [m/s]');
    legend('VelRaw','VertVel','VelRawSurf','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' poly fit ',num2str(polyTimePeriod),' s'],'interpreter','none');
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
    grid on
    grid minor
       
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f5, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_vertVelsAndPoly',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
    
    %% Plot average level velocities
   
    % Prepare for plotting
    
    scaleFact=0.5;
    
    addVals=0.25:0.5:max(altEdges);
    
    velRawMadj=(velRawMeans-nanmean(velRawMeans,2))*scaleFact+addVals';
    velMadj=(velMeans-nanmean(velMeans,2))*scaleFact+addVals';
    velUnlagMadj=(velUnlagMeans-nanmean(velUnlagMeans,2))*scaleFact+addVals';
    velCorrPMadj=(velCorrPMeans-nanmean(velCorrPMeans,2))*scaleFact+addVals';
    
    figure
    hold on
    scatter(velRawMeans(middleLevel,:),vertVelUnlag);
    xlabel('velRaw');
    ylabel('VertVelUnlag');
          
    f3=figure('DefaultAxesFontSize',12);
    set(f3,'Position',[200 500 3000 1500]);
    
    middleLevel=round(maxEdge/2);
    
    subplot(3,1,1)
    hold on
    plot(data.time,velRawMeans(middleLevel,:)-1,'-m');
    plot(data.time,velMeans(middleLevel,:),'-b');
    plot(data.time,velCorrPMeans(middleLevel,:)+1,'-k');
    plot(data.time,velCorrP30Means(middleLevel,:)+1,'-','color',[0.5 0 0]);
    plot(data.time,velCorrP60Means(middleLevel,:)+1,'-','color',[0 0.5 0]);
    plot(data.time,velCorrP120Means(middleLevel,:)+1,'-','color',[0 0 0.5]);
    plot(data.time,data.vertVel,'-c');
    plot(data.time,vertVelUnlag,'-g');
    plot(data.time,movmean(velRawS,50),'-r');
    plot(data.time,velRawMeans(middleLevel,:)-vertVelUnlag-0.1,'-','color',[0.5 0.5 0.5]);
    ylabel('Velocity [m/s]');
    ylim([-1 4]);
    
    yyaxis right
    plot(data.time,data.alt./1000,'-','color',[0.5 0 0]);
    ylabel('Altitude [km]');
    altLim=max(data.alt./1000)-min(data.alt./1000);
    ylim([min(data.alt./1000)-altLim*0.05 max(data.alt./1000)+altLim*2]);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
    grid on
    grid minor
    legend('VelRaw-1','Vel','VelCorrPoly+1','VelCorrPoly30+1','VelCorrPoly60+1','VelCorrPoly120+1',...
        'VertVel','VertVelUnlag','VelRawSurf','VelRaw-VertVelUnlag','Alt','Orientation','horizontal');
   
    title(['Vertical velocity around ',num2str(middleLevel),' km']);
    
    subplot(3,1,2)
    hold on
    plot(data.time,movmean(velRawS,50)-1,'-r');
    plot(data.time,movmean(velS,50),'-b');
    plot(data.time,movmean(velUnlagS,50),'-','color',[0.5 0.5 0.5]);
    plot(data.time,movmean(velCorrPS,50)+1,'-k');
    %ylim([-0.2 6]);
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
    grid on
    grid minor
    legend('VelRawSurf-1','VelSurf','VelUnlagSurf','VelCorrPoly3Surf+1','Orientation','horizontal');
   
    ylabel('Velocity [m/s]');
    title(['Vertical velocity of the ocean surface']);
    
    subplot(3,1,3)
    hold on
    plot(data.time,data.pitch,'-g');
    plot(data.time,data.tilt*-1,'-k');
    ylim([0 2]);
    ylabel('Pitch/tilt [deg]');
            
    yyaxis right
    plot(data.time,data.roll,'-b');
    plot(data.time,180-data.rot,'-c');
    %ylim([-0.2 6]);
    ylabel('Roll/Rot [deg]');
    xlim([startTime,endTime]);
    ax=gca;
    ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
    grid on
    grid minor
    legend('Pitch','-Tilt','Roll','180-Rot','Orientation','horizontal');
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f3, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_line'],'-dpng','-r0');
    end
    
    
   
    %% Plot vel field
      
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 1500 1500]);
    
    % Vel Raw
    subplot(4,1,1)
    hold on
    fig2=surf(data.time,asl,data.velRaw);
    fig2.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig2.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
    %limits=[-8 -2.25 -1.95 -1.65 -1.35 -1.05 -0.75 -0.45 -0.15 0.15 0.45 0.75 1.05 1.35 1.65 1.95 2.25 8];
    limits=[-8 -3.8 -3.3 -2.8 -2.3 -1.8 -1.3 -0.8 -0.3 0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 8];
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
    
    plot(data.time,velRawMadj,'-k');
    title('VelRaw');
    
    % Vel
    
    subplot(4,1,2)
    hold on
    fig2=surf(data.time,asl,data.vel);
    fig2.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig2.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
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
    
    plot(data.time,velMadj,'-k');
    
    title('Vel');
      
    % Vel  unlag
    
    subplot(4,1,3)
    hold on
    fig4=surf(data.time,asl,data.velRaw-vertVelUnlag);
    fig4.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig4.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
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
    
    plot(data.time,velUnlagMadj,'-k');
    
    title(['VelUnlag']);
    
    % Vel poly
    
    subplot(4,1,4)
    hold on
    fig4=surf(data.time,asl,velCorrP);
    fig4.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig4.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
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
    
    plot(data.time,velCorrPMadj,'-k');
       
    title(['VelCorr poly fit 3, ',num2str(polyTimePeriod),' s']);
       
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_vel'],'-dpng','-r0');
    end
    
    %% Test ideas
    %figure
    %scatter(velMeans(middleLevel,:),data.vertVel);
    
end