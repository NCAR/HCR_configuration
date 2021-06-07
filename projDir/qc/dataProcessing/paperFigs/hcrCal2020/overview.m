% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

loadData=1;

wi=10;
hi=8;

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
conc=imread('/h/eol/romatsch/papers/HCRcalibration/figs/HCR_concept.png');

s1=subplot(3,1,1);
imshow(conc);

%% Foto
foto=imread('/h/eol/romatsch/papers/HCRcalibration/figs/HCR_image5.jpg');
fotoSmall=foto(170:end,100:end,:);

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

s1.Position = [-0.01 0.59 1.03 0.6];
s2.Position=[0.01 0.46 0.4 0.36];
s3.Position=[0.475 0.445 0.5 0.31];
s4.Position=[0.06 0.065 0.36 0.38];
s5.Position=[0.475 0.065 0.5 0.31];

text(30,-450,'(a)','fontsize',12,'fontweight','bold');
text(30,52,'(b)','fontsize',12,'fontweight','bold','BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
text(360,450,'HCR','fontsize',12,'fontweight','bold','BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
annotation('arrow',[0.185 0.208],[0.585,0.571],'Linewidth',2,'color',[1 1 1]);

print(fig1, [figdir,'overview.png'],'-dpng','-r0');