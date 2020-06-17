% Compare data from before and after velocity correction

clear all
close all

savefig=0;

startTime=datetime(2018,2,7,21,30,0);
endTime=datetime(2018,2,7,21,35,0);
%endTime=datetime(2018,1,16,0,21,0);

addpath('/h/eol/romatsch/gitPriv/utils/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velDiffFigs_compare_polyfit/';
formatOut = 'yyyymmdd_HHMM';

fileList=makeFileList(indir,startTime,endTime,1);

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
        velCorrIn=ncread(infile,'VEL_CORR');
        data.velCorr=[data.velCorr,velCorrIn];
        data.velDiff1=[data.velDiff1,velCorrIn-velIn];
    end
    
    velDiff11d=nanmax(data.velDiff1);
    velDiff11d(find(velDiff11d==0))=nan;
   
    
    %% Get Scott's data
    [velCorrS,velCorrSurfS,rangeToSurfS,velSmoothS]=vel2vel_corr(data.vel,data.dbz,data.range,data.alt);
    
    %% Make poly fit
    polyTimePeriod=10; %Time period for poly fit in seconds
    polyOrder=3; % Order of polynomial fit
    [velCorrP,velSmoothP]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
     polyOrder=5; % Order of polynomial fit
    [velCorrP5,velSmoothP5]=vel2vel_corr_poly(data,polyTimePeriod,polyOrder);
    
    % Plot line plots
    
    close all
    
    f1=figure('DefaultAxesFontSize',14);
    set(f1,'Position',[1700 500 1500 500]);
    
    hold on
    plot(data.time,-velDiff11d,'-r','linewidth',1.2);
    plot(data.time,velSmoothS,'-k','linewidth',1.2);
    plot(data.time,velSmoothP,'-b','linewidth',1.2);
    plot(data.time,velSmoothP5,'-c','linewidth',1.2);
    ylim([-2 2]);
    ylabel('Vel smooth [m/s]');
    legend('VelSmooth CSET','VelSmooth Scott','VelSmooth Poly3','VelSmooth Poly5','Orientation','horizontal');
    title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut),' poly fit ',num2str(polyTimePeriod),' s'],'interpreter','none');
    xlim([startTime,endTime]);
    grid on
   
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_poly',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
    
    % Plot vel field
    
    %calculate above sea level altitudes
    asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 1500 1500]);
    
    subplot(4,1,1)
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
    title('Vel');
    
    subplot(4,1,2)
    fig2=surf(data.time,asl,data.velCorr);
    fig2.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig2.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
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
    title('VelCorr CSET');
    
    
    subplot(4,1,3)
    fig3=surf(data.time,asl,velCorrS);
    fig3.EdgeColor='none';
    ylim([-0.2 6]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig3.CData;
    
    col_def1 = nan(size(fld));
    col_def2 = nan(size(fld));
    col_def3 = nan(size(fld));
    
    color_map=colormap(vel_default);
    
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
    
    subplot(4,1,4)
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
    title(['VelCorr poly fit 3, ',num2str(polyTimePeriod),' s']);
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_vel_',num2str(polyTimePeriod),'s'],'-dpng','-r0');
    end
    
end