% Analyze HCR clouds

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
addName=''; % Extra name part for output files. Default is ''.

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];
caseList = table2array(readtable(infile));

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

figdir=['/h/eol/romatsch/hcrCalib/mask/',project,'/'];
outdir=['/scr/sci/romatsch/data/HCRprods/mask/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

directories.dataDir=HCRdir(project,quality,freqData);

cloudInds=[];
timeFlight=[];

%% Go through flight hours
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
    
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
        
        % Desired variables. The variable name comies after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        %% Load data
        
        data.DBZ=[];
        %data.VEL=[];
        %data.VEL_RAW=[];
        %data.VEL_CORR=[];
        data.WIDTH=[];
        %data.WIDTH_CORR=[];
        %data.DBMVC=[];
        %data.SNR=[];
        %data.NCP=[];
        %data.LDR=[];
        
        dataVars=fieldnames(data);
        
        % Make list of files within the specified time frame
        fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            startTime=endTime;
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        timeFlight=cat(2,timeFlight,data.time);
        %% Mask data
        
        [maskData antStat]=echoMask(data);
        
        refl=data.DBZ;
        refl(maskData>1)=nan;
        
        maskPlot=maskData;
        maskPlot(maskData==0)=nan;
        
        cloudIndsHour=zeros(size(antStat));
        cloudIndsHour(any(maskData==1,1))=1;
        cloudIndsHour(any(maskData==8,1))=1;
        
        cloudInds=cat(2,cloudInds,cloudIndsHour);
        %% Plot
        ytickLabels={'Good echo';'Bang';'Surface';'Below surf.';...
            'Backlobe';'NS cal';'Ant. trans.';'Missing'};
        
        colMask=[0,0,1;
            1,0,0;
            0.5,0.5,0.5;
            0.7065,0.4415,0.2812;
            0,1,0;
            0.9290,0.8940,0.1250;
            1,0,1;
            0,0,0];
        
        close all
        
        figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200]);
        
        ax1=subplot(3,1,1);
        fig1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
        view(2);
        fig1=colMapDBZ(fig1);
        ylim([-1 8]);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        
        ax3=subplot(3,1,3);
        fig3=surf(data.time,data.asl./1000,refl,'edgecolor','none');
        view(2);
        fig3=colMapDBZ(fig3);
        ylim([-1 8]);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        
        ax2=subplot(3,1,2);
        fig2=surf(data.time,data.asl./1000,maskPlot,'edgecolor','none');
        view(2);
        caxis([1 8]);
        colormap(ax2,colMask);
        hcb=colorbar;
        set(hcb,'ytick',[1.5:0.86:8.5]);
        set(hcb,'YTickLabel',ytickLabels);
        ylim([-1 8]);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        
        %linkaxes([ax1,ax2,ax3],'xy');
        
        formatOut = 'yyyymmdd_HHMM';
        set(gcf,'PaperPositionMode','auto')
        print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_echoMask'],...
            '-dpng','-r0');
        
        startTime=endTime;
        
        time=data.time;
        save([outdir,project,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.mat'],'time');
        
        mask=maskData;
        save([outdir,project,'.mask.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.mat'],'mask');
        
        save([outdir,project,'.antStat.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.mat'],'antStat');
    end
    
    save([outdir,project,'.cloudInds.',datestr(timeFlight(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(timeFlight(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'cloudInds');
    save([outdir,project,'.timeFlight.',datestr(timeFlight(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(timeFlight(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'timeFlight');
end