% Sum up flight hours per grid cell

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

projects={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersions={'v3.0','v3.1','v3.1'};

gridStep=1;

savedir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/';
flightHourGrids={};

%% Set up grids

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

for bb=1:length(projects)
    project=projects{bb};
    qcVersion=qcVersions{bb};

    disp(project);

    lonLength=(lonLims(bb,2)-lonLims(bb,1))/gridStep;
    latLength=(latLims(bb,2)-latLims(bb,1))/gridStep;
    hourGrid=zeros(latLength,lonLength);

    lonSteps=lonLims(bb,1):gridStep:lonLims(bb,2);
    latSteps=latLims(bb,1):gridStep:latLims(bb,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cfDir=HCRdir(project,quality,qcVersion,freqData);

    infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

    caseList=table2array(readtable(infile));


    %% Loop through flights

    for aa=1:size(caseList,1)
        disp(['Flight ',num2str(aa)]);

        startTime=datetime(caseList(aa,1:6));
        endTime=datetime(caseList(aa,7:12));

        fileList=makeFileList(cfDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

        data=[];
        data.dummy=[];

        % Load data
        data=read_HCR(fileList,data,startTime,endTime);

        %% Loop through output grid
        for ii=1:size(hourGrid,1)
            for jj=1:size(hourGrid,2)
                pixInds=find(data.latitude>latSteps(ii) & data.latitude<=latSteps(ii+1) & ...
                    data.longitude>lonSteps(jj) & data.longitude<=lonSteps(jj+1));
                hzPix=length(pixInds);
                hourGrid(ii,jj)=hourGrid(ii,jj)+hzPix/10/60/60;
            end
        end

    end

    flightHourGrids.(project)=hourGrid;

end

%% Save properties

disp('Saving output ...');

save([savedir,'flightHourGrids.mat'],'flightHourGrids');
