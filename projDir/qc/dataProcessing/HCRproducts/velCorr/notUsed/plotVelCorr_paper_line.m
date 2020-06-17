% Compare data from before and after velocity correction

clear all
close all

fig=5; %1, 5, 8, 10, 12, 101
savefig=1;

addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/utils/colormaps/');

formatOut = 'yyyymmdd_HHMM';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velFigs_paper/';

if fig==1
    startTime=datetime(2015,2,2,19,15,0);
    endTime=datetime(2015,2,2,19,20,0);
    minEdge1=-3;
    maxEdge1=1.5;
    minEdge2=20;
    maxEdge2=80;
end

if fig==5
    startTime=datetime(2015,2,2,19,15,0);
    endTime=datetime(2015,2,2,19,20,0);
    minEdge1=-1.5;
    maxEdge1=1.5;
end

if fig==8
    startTime=datetime(2015,7,17,17,35,0);
    endTime=datetime(2015,7,17,17,45,0);
    indir='/scr/eldora2/rsfdata/cset/hcr/qc/cfradial/velcorr/10hz/';
    xtoff=15;
    ytoff=0.05;
    minEdge1=-1.5;
    maxEdge1=1.5;
end

if fig==10
    startTime=datetime(2018,1,22,23,5,0);
    endTime=datetime(2018,1,22,23,20,0);
    xtoff=15;
    ytoff=0.05;
    minEdge1=-1.5;
    maxEdge1=1.5;
end

if fig==12
    startTime=datetime(2018,1,23,1,50,0);
    endTime=datetime(2018,1,23,2,0,0);
    xtoff=15;
    ytoff=0.05;
    minEdge1=-1.5;
    maxEdge1=1.5;
end

if fig==101
    startTime=datetime(2015,1,31,18,58,0);
    endTime=datetime(2015,1,31,18,59,0);
    minEdge1=-3;
    maxEdge1=1.5;
    minEdge2=20;
    maxEdge2=80;
end

if fig==1 | fig==5
    indir='/scr/snow2/rsfdata/projects/noreaster/cfradial/moments/qcv1/10hz/';
end

if fig==101
    indir='/scr/eldora1/rsfdata/noreaster/cfradial/moments/orig/10hz/';
end

if fig==10 | fig==12
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
        
        if fig~=101
            velCorrIn=ncread(infile,'VEL_CORR');
            data.velCorr=[data.velCorr,velCorrIn];
        end
    end
    
    %% Get surface indices
    % Find surface indices
    DBZ=data.dbz;
    DBZ(1:15,:)=0;
    
    rangeTemp=data.range;
    DBZmask=DBZ; %Mask with DBZ values that are close to the surface
    
    rangeTemp(find(abs(data.range-data.alt)>2000))=nan;
    DBZmask(isnan(rangeTemp))=nan;
    [bla ground_index]=nanmax(DBZmask,[],1);
    
    % Convert to linear indices
    linInd=sub2ind(size(DBZ),ground_index,1:length(data.alt));
    
    surfDBZ=data.dbz(linInd);
    surfVel=data.vel(linInd);
    if fig~=101
        surfVelCorr=data.velCorr(linInd);
    end
    
    %% Plot
    
    if fig==1 | fig==101
        
        close all
        
        hm=2;
        wm=2;
        wi=10;
        hi=4;
        
        fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
        fig1.PaperPositionMode = 'manual';
        fig1.PaperUnits = 'inches';
        fig1.Units = 'inches';
        fig1.PaperPosition = [0, 0, wi, hi];
        fig1.PaperSize = [wi, hi];
        fig1.Resize = 'off';
        fig1.InvertHardcopy = 'off';
        
        set(fig1,'color','w');
        
        ax1=subtightplot(1,1,1,[hm,wm]);
        hold on;
        outerpos1 = ax1.Position;
        ax1.Position = [outerpos1(1)+0.03 outerpos1(2)+0.07 outerpos1(3)-0.04 outerpos1(4)-0.04];
        
        hold on;
        plot(data.time,surfVel,'-b','linewidth',1.2);
        xlim([startTime,endTime]);
        ylim([minEdge1 maxEdge1]);
        ylabel('v_{surf}^{meas} [m s^{-1}]');
        
        yyaxis right
        plot(data.time,surfDBZ,'-r','linewidth',1.2);
        xlim([startTime,endTime]);
        ylim([minEdge2 maxEdge2]);
        ylabel('Z_e [dBZ]');
        
        plt = gca;
        plt.YAxis(2).Color = 'k'; % change color of RHS y-axis to black
        
        legend('v_{surf}^{meas}','Z_e');
        %text(startTime+seconds(xtoff),maxEdge1-ytoff,'a','fontsize',16,'fontweight','bold');
    end
    
    if fig==5
        
        close all
       
        hm=2;
        wm=2;
        wi=10;
        hi=3;
        
        fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
        fig1.PaperPositionMode = 'manual';
        fig1.PaperUnits = 'inches';
        fig1.Units = 'inches';
        fig1.PaperPosition = [0, 0, wi, hi];
        fig1.PaperSize = [wi, hi];
        fig1.Resize = 'off';
        fig1.InvertHardcopy = 'off';
        
        set(fig1,'color','w');
        
        ax1=subtightplot(1,1,1,[hm,wm]);
        hold on;
        outerpos1 = ax1.Position;
        ax1.Position = [outerpos1(1)+0.03 outerpos1(2)+0.1 outerpos1(3)-0.02 outerpos1(4)-0.08];
        
        hold on;
        plot(data.time,surfVelCorr,'-k','linewidth',1.2);
        xlim([startTime,endTime]);
        ylim([minEdge1 maxEdge1]);
        ylabel('v_{surf}^{corr} [m s^{-1}]');
        
        text(datetime(2015,2,2,19,15,25),1.17,[{'Mean: -0.004 m s^{-1}'};{'Variance: 0.008 m^2 s^{-2}'}],'fontsize',12);
               
        legend('v_{surf}^{corr}');
        %text(startTime+seconds(xtoff),maxEdge1-ytoff,'a','fontsize',16,'fontweight','bold');
    end
    
    if fig==8 | fig==10 | fig==12
        
        close all
       
        hm=0.2;
        wm=0.2;
        wi=10;
        hi=5;
        
        fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
        fig1.PaperPositionMode = 'manual';
        fig1.PaperUnits = 'inches';
        fig1.Units = 'inches';
        fig1.PaperPosition = [0, 0, wi, hi];
        fig1.PaperSize = [wi, hi];
        fig1.Resize = 'off';
        fig1.InvertHardcopy = 'off';
        
        set(fig1,'color','w');
        
        ax1=subtightplot(2,1,1,[hm,wm]);
        hold on;
        outerpos1 = ax1.Position;
        ax1.Position = [outerpos1(1)+0.03 outerpos1(2)+0.0 outerpos1(3)-0.02 outerpos1(4)+0.02];
        
        hold on;
        leg1=plot(data.time,surfVel,'-b','linewidth',0.8);
        leg2=plot(data.time,movmean(surfVel,50),'-r','linewidth',1.5);
        xlim([startTime,endTime]);
        ylim([minEdge1 maxEdge1]);
        ylabel('v_{surf}^{meas} [m s^{-1}]');
               
        text(startTime+seconds(xtoff),maxEdge1-ytoff,'a','fontsize',16,'fontweight','bold');
        
        if fig==12
            plot([datetime(2018,1,23,1,56,30) datetime(2018,1,23,1,56,30)],[minEdge1, maxEdge1],'-g','linewidth',2);
            plot([datetime(2018,1,23,1,57,25) datetime(2018,1,23,1,57,25)],[minEdge1, maxEdge1],'-g','linewidth',2);
            ax = gca;
            ax.SortMethod = 'childorder';
        end
        legend([leg1,leg2],{'v_{surf}^{meas}','5s run mean'});
        
        ax1=subtightplot(2,1,2,[hm,wm]);
        hold on;
        outerpos1 = ax1.Position;
        ax1.Position = [outerpos1(1)+0.03 outerpos1(2)+0.05 outerpos1(3)-0.02 outerpos1(4)+0.02];
        
        hold on;
        leg1=plot(data.time,surfVelCorr,'-b','linewidth',0.8);
        leg2=plot(data.time,movmean(surfVelCorr,50),'-r','linewidth',1.5);
        xlim([startTime,endTime]);
        ylim([minEdge1 maxEdge1]);
        ylabel('v_{surf}^{corr} [m s^{-1}]');
        
        text(startTime+seconds(xtoff),maxEdge1-ytoff,'b','fontsize',16,'fontweight','bold');
        
        if fig==12
            plot([datetime(2018,1,23,1,56,30) datetime(2018,1,23,1,56,30)],[minEdge1, maxEdge1],'-g','linewidth',2);
            plot([datetime(2018,1,23,1,57,25) datetime(2018,1,23,1,57,25)],[minEdge1, maxEdge1],'-g','linewidth',2);
            ax = gca;
            ax.SortMethod = 'childorder';
        end
        legend([leg1,leg2],{'v_{surf}^{corr}','5s run mean'});
    end
    
    if savefig
        set(gcf,'PaperPositionMode','auto')
        print(fig1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_fig',num2str(fig)],'-dpng','-r0');
    end
    
end