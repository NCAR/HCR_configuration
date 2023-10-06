% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates (Frozen and Liquid), cset (Liquid)
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

phase='Liquid';

showPlot='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

%figdir=[dataDir(1:end-5),'attenPlots/findCoeffs/'];
figdir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.2_full/attenPlots/findCoeffs/';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/atten',phase,'_',project,'_dev.txt'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    %% Get data

    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ = [];
    data.U=[];
    data.V=[];
    data.SST=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    data.ANTFLAG=[];
    data.rotation=[];
    data.MELTING_LAYER=[];
    data.pulse_width=[];

    dataVars=fieldnames(data);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.frq=ncread(fileList{1},'frequency');

    %% Get surface wind
    [linInd,~,~]=hcrSurfInds(data);

    Utemp=fillmissing(data.U,'previous',1);
    Vtemp=fillmissing(data.V,'previous',1);
    data.U_SURF=Utemp(linInd);
    data.V_SURF=Vtemp(linInd);

    %% Correct for gaseous attenuation

    disp('Calculating gaseous attenuation ...');
    [~,gasAttCloud,~,gasAttCloudMat]=get_gas_atten(data);

    % Extend to surface
    
    for jj=1:2
        replaceLinInd=linInd;
        replaceLinInd(isnan(linInd) | ~isnan(gasAttCloudMat(linInd)) | data.FLAG(linInd)~=7)=[];

        gasAttCloudMat(replaceLinInd)=gasAttCloudMat(replaceLinInd-jj);
        if jj==2
            gasAttCloudMat(replaceLinInd-1)=gasAttCloudMat(replaceLinInd-jj);
        else
            replaceLinInd1=replaceLinInd;
        end
    end

    [~,cInd]=ind2sub(size(data.DBZ),replaceLinInd1);
    gasAttCloud(cInd)=gasAttCloud(cInd)'+gasAttCloudMat(replaceLinInd1);
    
    %PIA gas
    piaGasMat2=cumsum(gasAttCloudMat,1).*2;

    data.DBZcorrGas=data.DBZ+piaGasMat2;

    %% Remove all up pointing and unsuitable data

    dbzOrig=data.DBZ;

    % Noise source cal (10), missing (11)
    badInds=find(any(data.FLAG>9,1));

    % Land surface
    badInds=cat(2,badInds,find(any(data.FLAG==8,1)));
    
    % Zenith (2), pointing (3), scanning (4), transision (5), failure (6)
    badInds=cat(2,badInds,find(data.ANTFLAG>1));

    infields=fields(data);
    for bb=1:length(infields)
        if strcmp(infields{bb},'DBZ') | strcmp(infields{bb},'FLAG') | ...
                strcmp(infields{bb},'rotation') | strcmp(infields{bb},'elevation')
            currfield=data.(infields{bb});
            currfield(:,badInds)=nan;
            data.(infields{bb})=currfield;
        end
    end

    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG~=1)=nan;

    data.dbzMaskedCorrGas=data.DBZcorrGas;
    data.dbzMaskedCorrGas(data.FLAG~=1)=nan;

   %% Get hydrometeor attenuation

   [sig0measAtt,surfFlag,refSig0,sig0model,piaHydromet2]=getRefAtten_fromGasCorr(data);
   
   %% Use only data of clouds with precip to the ground and no big gaps
    testDBZ=data.dbzMaskedCorrGas;
    
    testDBZ(1:30,:)=nan;

    firstInd=nan(1,size(testDBZ,2));
    for ii=1:size(testDBZ,2)
        thisRay=testDBZ(:,ii);
        if any(~isnan(thisRay))
            firstInd(ii)=find(~isnan(thisRay),1,'first');
            lastInd=find(~isnan(thisRay),1,'last');
            firstSurf=find(data.FLAG(:,ii)==7,1,'first');
            if lastInd+1~=firstSurf
                testDBZ(:,ii)=nan;
                continue
            else
                cloudRay=thisRay(firstInd(ii):lastInd);
                nanNum=length(find(isnan(cloudRay)));
                if nanNum/length(cloudRay)>0.05
                    testDBZ(:,ii)=nan;
                end
            end
        end
    end

   %% Hitschfeld Bordan from surface up
   piaHydromet1=piaHydromet2/2;

   % Fix PIA
   piaHydrometInt=piaHydromet1;
   piaHydrometInt(piaHydrometInt<0)=0;

   % Handle extinct
   extInds=any(data.FLAG==3,1);
   extInds=imdilate(extInds,strel('disk',3));
   piaHydrometInt(isnan(piaHydrometInt) & extInds==0)=0;

   piaHydrometInt=fillmissing(piaHydrometInt,'linear','EndValues','nearest');

   cloudInds=any(data.FLAG==1,1);
   piaHydrometInt(~cloudInds)=nan;
   piaHydrometInt(badInds)=nan;

   if strcmp(phase,'Frozen')
    abGuess=[0.0325,1];
   elseif strcmp(phase,'Liquid')
       abGuess=[0.05,1];
   end

   %% Optimize
   disp('Optimizing ...');

   fun=@(abGuess) hitschfeldBordan_optimize(testDBZ,piaHydrometInt,data.range,abGuess,firstInd);
   [abOpt,fval,exitflag,output] = fminsearch(fun,abGuess);

   disp(['alpha=',num2str(abOpt(1)),', beta=',num2str(abOpt(2))]);

   zHB=hitschfeldBordan_newCoeffs(testDBZ,piaHydrometInt,data.range,abOpt);

   attCorrection=zHB-data.dbzMaskedCorrGas;

   numNeg=length(find(attCorrection<0));
   minNeg=min(attCorrection(:),[],'omitnan');

   disp(['Negative pixels: ',num2str(numNeg),', minimum: ',num2str(minNeg),' dB']);

   goodInds=find(~isnan(firstInd));
   topZdiff=[];
   for ii=1:length(goodInds)
       testZ=attCorrection(firstInd(goodInds(ii)),ii);
       if ~isnan(testZ)
           topZdiff=[topZdiff,attCorrection(firstInd(goodInds(ii)),ii)];
       end
   end

   ftest = figure('Position',[200 500 1800 300],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);
   plot(topZdiff);

 
   %% Plot

   timeInds=1:round(length(data.time)/2000):length(data.time);

   ylimRefl=ceil(max(data.asl(~isnan(zHB)))./1000);

   sig0measClear=nan(size(data.time));
   sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
   sig0measCloud=nan(size(data.time));
   sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);

   refSig0(badInds)=nan;
 
   close all

   f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);

   colormap dbz_default

   s1=subplot(5,1,1);
   hold on
   surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,dbzOrig(:,timeInds),'edgecolor','none');
   view(2);
   ylabel('Altitude (km)');
   caxis([-20 30]);
   ylim([-0.1 ylimRefl]);
   xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
   colorbar
   grid on
   title('Reflectivity (dBZ)')

   s2=subplot(5,1,2);
   hold on
   surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.dbzMaskedCorrGas(:,timeInds),'edgecolor','none');
   view(2);
   ylabel('Altitude (km)');
   caxis([-20 30]);
   ylim([-0.1 ylimRefl]);
   xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
   colorbar
   grid on
   title('Reflectivity corrected for gaseous attenuation (dBZ)')

   s3=subplot(5,1,3);
   hold on
   surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,zHB(:,timeInds),'edgecolor','none');
   view(2);
   ylabel('Altitude (km)');
   caxis([-20 30]);
   ylim([-0.1 ylimRefl]);
   xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
   colorbar
   grid on
   title('Reflectivity corrected for gaseous and liquid attenuation (dBZ)')

   s4=subplot(5,1,4);
   colmapCorr=cat(1,[1,0,0],turbo(15));
   hold on
   surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,attCorrection(:,timeInds),'edgecolor','none');
   view(2);
   ylabel('Altitude (km)');
   caxis([-1.01 15]);
   s4.Colormap=colmapCorr;
   ylim([-0.1 ylimRefl]);
   xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
   colorbar
   grid on
   title(['Attenuation correction (dBZ), alpha=',num2str(abOpt(1)),', beta=',num2str(abOpt(2))]);

   s5=subplot(5,1,5);
   hold on
   l0=plot(data.time,sig0model,'-c','linewidth',2);
   l1=plot(data.time,sig0measClear,'-b','linewidth',1);
   l2=plot(data.time,sig0measCloud,'color',[0.5 0.5 0.5],'linewidth',0.5);
   l3=plot(data.time,refSig0,'-r','linewidth',2);
   ylabel('Sig0 (dB)');
   ylim([0 20]);

   yyaxis right
   l6=plot(data.time,gasAttCloud*2,'-k','linewidth',1);
   l7=plot(data.time,piaHydrometInt*2,'-','color',[0,0.5,0],'linewidth',2);
   l8=plot(data.time,piaHydromet2,'-g','linewidth',1);
   ylabel('Atten. (dB)');
   ylim([-5 15]);
   grid on
   set(gca,'YColor','k');

   xlim([data.time(1),data.time(end)]);

   legend([l0 l1 l2 l3 l6 l7 l8],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
       'sig0 ref','2-way gas att','2-way liq att int','2-way liq att'},...
       'orientation','horizontal','location','north');
   title([datestr(data.time(1),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(end),'yyyy-mm-dd HH:MM:SS')])

   s1Pos=s1.Position;
   s2Pos=s2.Position;
   s3Pos=s3.Position;
   s4Pos=s4.Position;
   s5Pos=s5.Position;

   s1.Position=[s1Pos(1),s1Pos(2),s5Pos(3),s1Pos(4)];
   s2.Position=[s2Pos(1),s2Pos(2),s5Pos(3),s2Pos(4)];
   s3.Position=[s3Pos(1),s3Pos(2),s5Pos(3),s3Pos(4)];
   s4.Position=[s4Pos(1),s4Pos(2),s5Pos(3),s4Pos(4)];

   set(gcf,'PaperPositionMode','auto')
   print(f1,[figdir,project,'_att',phase,'_',...
       datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')

end