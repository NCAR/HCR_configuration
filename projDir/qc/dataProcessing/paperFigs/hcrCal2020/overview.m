% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

loadData=1;

wi=10;
hi=10;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%% Conceptual drawing
conc=imread('/h/eol/romatsch/papers/HCRcalibration/figs/HCR.jpg');
concSmall=conc(900:1330,700:2850,:);

s1=subplot(3,1,1);
imshow(concSmall);

%% Foto
foto=imread('/h/eol/romatsch/papers/HCRcalibration/figs/HCR_image5.jpg');
fotoSmall=foto(170:end,300:end,:);

s2=subplot(3,2,3);
imshow(fotoSmall);

%% CSET

s3=subplot(3,2,4);
hold on;

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
qcVersion='v2.1';
freq='10hz';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];
caseList = table2array(readtable(infile));
indir=HCRdir(project,quality,qcVersion,freq);

xlims=[-165 -110];
ylims=[15 45];

xlim(xlims);
ylim(ylims);

title('(c) CSET');
%xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

S = shaperead('landareas','UseGeoCoords',true);
geoshow([S.Lat],[S.Lon],'Color','k')
%geoshow(S ,'FaceColor', [0.8 0.8 0.8])
box on

if loadData
    for ii=1:size(caseList,1)
        disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
        
        startTime=datetime(caseList(ii,1:6));
        endTime=datetime(caseList(ii,7:12));
        
        data=[];
        data.dummy=[];
        
        %% Load data
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            startTime=endTime;
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        %% Plot
        
        plot(data.longitude,data.latitude,'-b');
    end    
end

%% SOCRATES

s4=subplot(3,2,5);
hold on;

project='socrates'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
qcVersion='v2.1';
freq='10hz';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];
caseList = table2array(readtable(infile));
indir=HCRdir(project,quality,qcVersion,freq);

xlims=[120 180];
ylims=[-70 -30];

xlim(xlims);
ylim(ylims);

title('(d) SOCRATES');
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

geoshow([S.Lat],[S.Lon],'Color','k')
box on

if loadData
    for ii=1:size(caseList,1)
        disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
        
        startTime=datetime(caseList(ii,1:6));
        endTime=datetime(caseList(ii,7:12));
        
        data=[];
        data.dummy=[];
        
        %% Load data
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            startTime=endTime;
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        %% Plot
        
        plot(data.longitude,data.latitude,'-b');
    end    
end

%% OTREC

s5=subplot(3,2,6);
hold on;

project='otrec'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
qcVersion='v2.2';
freq='10hz';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];
caseList = table2array(readtable(infile));
indir=HCRdir(project,quality,qcVersion,freq);

xlims=[-100 -70];
ylims=[0 15];

xlim(xlims);
ylim(ylims);

title('(e) OTREC');
xlabel('Longitude (deg)')
%ylabel('Latitude (deg)')

geoshow([S.Lat],[S.Lon],'Color','k')

box on

if loadData
    for ii=1:size(caseList,1)
        disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
        
        startTime=datetime(caseList(ii,1:6));
        endTime=datetime(caseList(ii,7:12));
        
        data=[];
        data.dummy=[];
        
        %% Load data
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            startTime=endTime;
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        %% Plot
        data.latitude(data.longitude>-75)=[];
        data.longitude(data.longitude>-75)=[];
        
        plot(data.longitude,data.latitude,'-b');
    end    
end

s1.Position = [0.01 0.59 0.98 0.6];
s2.Position=[0.01 0.43 0.4 0.36];
s3.Position=[0.475 0.42 0.5 0.31];
s4.Position=[0.06 0.05 0.36 0.38];
s5.Position=[0.475 0.05 0.5 0.31];

text(30,-550,'(a)','fontsize',12,'fontweight','bold');
text(30,52,'(b)','fontsize',12,'fontweight','bold','BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(200,450,'HCR','fontsize',12,'fontweight','bold','BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
annotation('arrow',[0.145 0.172],[0.555,0.545],'Linewidth',2,'color',[1 1 1]);

text(10,-60,'Cooling exhaust','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(310,-60,'Data system','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(560,-60,'Klystron power','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(900,-60,'Radar front end','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(1200,-60,'Cooling inlet','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(1460,-60,'Antenna','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(1750,-60,'Reflector','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(2000,-60,'Inertial nav head','fontsize',10,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);

print(fig1, [figdir,'overview.png'],'-dpng','-r0');