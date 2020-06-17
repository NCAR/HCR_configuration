% Compare data from before and after velocity correction

clear all
close all

fig=9; %4, 6, 7, 9, 11
savefig=1;

addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/utils/colormaps/');

formatOut = 'yyyymmdd_HHMM';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velFigs_paper/';

if fig==4
    startTime=datetime(2015,2,2,19,15,0);
    endTime=datetime(2015,2,2,19,20,0);
    xtoff=7;
    ytoff=0.6;
    maxEdge=9;
end

if fig==6
    startTime=datetime(2015,2,2,14,5,0);
    endTime=datetime(2015,2,2,14,10,0);
    xtoff=7;
    ytoff=0.5;
    maxEdge=10;
end

if fig==7
    startTime=datetime(2015,7,17,17,35,0);
    endTime=datetime(2015,7,17,17,45,0);
    indir='/scr/eldora2/rsfdata/cset/hcr/qc/cfradial/velcorr/10hz/';
    xtoff=10;
    ytoff=0.05;
    maxEdge=1.2;
end

if fig==9
    startTime=datetime(2018,1,22,23,5,0);
    endTime=datetime(2018,1,22,23,20,0);
    xtoff=12;
    ytoff=1;
    maxEdge=5.5;
end

if fig==11
    startTime=datetime(2018,1,23,1,50,0);
    endTime=datetime(2018,1,23,2,0,0);
    xtoff=7;
    ytoff=0.25;
    maxEdge=3.5;
end

if fig==4 | fig==6
    indir='/scr/eldora1/rsfdata/noreaster/cfradial/moments/qcv1/10hz/';
end

if fig==9 | fig==11
    indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
end

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if ~isempty(fileList)
    
    data.time=[];
    data.alt=[];
    data.elev=[];
    data.range=[];
    data.dbz=[];
    data.vel=[];
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
        
        data.dbz=[data.dbz,ncread(infile,'DBZ')];
        
        velIn=ncread(infile,'VEL');
        data.vel=[data.vel,velIn];
        velCorrIn=ncread(infile,'VEL_CORR');
        data.velCorr=[data.velCorr,velCorrIn];
    end
    %% Remove bang
    if fig==9 | fig==11
        bangInds=find(data.range<70);
        data.dbz(bangInds)=nan;
        data.vel(bangInds)=nan;
        data.velCorr(bangInds)=nan;
    end
    
    %% Plot vel field
    close all
    %calculate above sea level altitudes
    asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
        
    hm=0.08;
    wm=0.1;
    wi=10;
    hi=10;
    
    fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
    fig1.PaperPositionMode = 'manual';
    fig1.PaperUnits = 'inches';
    fig1.Units = 'inches';
    fig1.PaperPosition = [0, 0, wi, hi];
    fig1.PaperSize = [wi, hi];
    fig1.Resize = 'off';
    fig1.InvertHardcopy = 'off';
    
    set(fig1,'color','w');
    
    colLines=lines;
    colClass=([0 1 1;colLines(2:4,:)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subtightplot(3,1,1,[hm,wm]);
    hold on;
    outerpos1 = ax1.Position;
    ax1.Position = [outerpos1(1)+0.02 outerpos1(2) outerpos1(3) outerpos1(4)];
    fig2=surf(data.time,asl,data.dbz);
    fig2.EdgeColor='none';
    ylim([-0.2 maxEdge]);
    xlim([startTime,endTime]);
    view(2);
    
    fld=fig2.CData;
    
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
    title('Z_e');
    text(startTime+seconds(xtoff),maxEdge-ytoff,'a','fontsize',16,'fontweight','bold');
    if fig==11
        plot([datetime(2018,1,23,1,56,30) datetime(2018,1,23,1,56,30)],[-0.2, maxEdge],'-g','linewidth',2);
        plot([datetime(2018,1,23,1,57,25) datetime(2018,1,23,1,57,25)],[-0.2, maxEdge],'-g','linewidth',2);
        ax = gca; 
        ax.SortMethod = 'childorder';
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subtightplot(3,1,2,[hm,wm]);
    hold on;
    outerpos1 = ax1.Position;
    ax1.Position = [outerpos1(1)+0.02 outerpos1(2) outerpos1(3) outerpos1(4)];
    fig2=surf(data.time,asl,data.vel);
    fig2.EdgeColor='none';
    ylim([-0.2 maxEdge]);
    xlim([startTime,endTime]);
    view(2);
    
    limits=[-8 -3.8 -3.3 -2.8 -2.3 -1.8 -1.3 -0.8 -0.3 0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 8];
    color_map=colormap(vel_default(17));
    
    %limits=-4.05:0.3:4.05;
    %color_map=colormap(vel_default(29));
    
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
    set(get(hcb,'Title'),'String','m s^{-1}');
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
    title('Vr^{meas}');
    text(startTime+seconds(xtoff),maxEdge-ytoff,'b','fontsize',16,'fontweight','bold');
    if fig==11
        plot([datetime(2018,1,23,1,56,30) datetime(2018,1,23,1,56,30)],[-0.2, maxEdge],'-g','linewidth',2);
        plot([datetime(2018,1,23,1,57,25) datetime(2018,1,23,1,57,25)],[-0.2, maxEdge],'-g','linewidth',2);
        ax = gca; 
        ax.SortMethod = 'childorder';
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VELcorr %%%%%%%%%%%%%%%%%%%%%%%
    ax1=subtightplot(3,1,3,[hm,wm]);
    hold on;
    outerpos1 = ax1.Position;
    ax1.Position = [outerpos1(1)+0.02 outerpos1(2) outerpos1(3) outerpos1(4)];
    fig3=surf(data.time,asl,data.velCorr);
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
    set(get(hcb,'Title'),'String','m s^{-1}');
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
    title('Vr^{corr}');
    text(startTime+seconds(xtoff),maxEdge-ytoff,'c','fontsize',16,'fontweight','bold');
    if fig==11
        plot([datetime(2018,1,23,1,56,30) datetime(2018,1,23,1,56,30)],[-0.2, maxEdge],'-g','linewidth',2);
        plot([datetime(2018,1,23,1,57,25) datetime(2018,1,23,1,57,25)],[-0.2, maxEdge],'-g','linewidth',2);
        ax = gca; 
        ax.SortMethod = 'childorder';
    end
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(fig1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_fig',num2str(fig)],'-dpng','-r0');
    end
    
end