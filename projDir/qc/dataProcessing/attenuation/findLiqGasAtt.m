% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

showPlot='on';
ylimRefl=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'attenPlots/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/atten_',project,'.txt'];

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
    data.U_SURF=[];
    data.V_SURF=[];
    data.SST=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    data.rotation=[];
    data.MELTING_LAYER=[];
    data.pulse_width=[];
            
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));

    data.frq=ncread(fileList{1},'frequency');

     %% Remove all up pointing data
    
    upInds=find(data.elevation>-85);
    upInds=cat(2,upInds,find(any(data.FLAG>9,1)),find(any(data.FLAG==3,1)));
    
    infields=fields(data);
    for bb=1:length(infields)
        if strcmp(infields{bb},'DBZ') | strcmp(infields{bb},'FLAG') | ...
                strcmp(infields{bb},'rotation') | strcmp(infields{bb},'elevation')
            currfield=data.(infields{bb});
            currfield(:,upInds)=nan;
            data.(infields{bb})=currfield;
        end
    end
    
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;

    %% One way and two way gaseous attenuation
    
    disp('Calculating gaseous attenuation ...');

    [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_gas_atten(data);
    data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY=2*gasAttCloud';

    %% Calculate sigma0 from model and from reflectivity
    
    disp('Calculating sig0 ...');
    
    % Find ocean surface gate
    [linInd,maxGate,rangeToSurf] = hcrSurfInds(data);
    
    % Measured sig0 from surface reflectivity
    data.surfRefl=data.DBZ(linInd);
    sig0measured=calc_sig0_surfRefl(data);
    
    sig0measAtt=sig0measured(linInd)+data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY;
    sig0measAtt(data.elevation>-85)=nan;
    
%     sig0measLin=10.^(sig0measured./10);
%     sig0meas3gates=nan(size(data.time));
%     for kk=1:length(data.time)
%         if ~isnan(maxGate(kk))
%             sig0meas3gates(kk)=sum(sig0measLin(maxGate(kk)-1:maxGate(kk)+1,kk),'omitnan');
%         end
%     end
%     
%     clear sig0measLin sig0measured
%     
%     sig0measAtt3gates=10.*log10(sig0meas3gates)+data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY;  
%     sig0measAtt3gates(data.elevation>0)=nan;
%     sig0measAtt3gates(imag(sig0measAtt3gates)~=0)=nan;
    
    % sig0 from models
    sig0modelAll= calc_sig0_model(data);
    %sig0modelFV=sig0modelAll(2,:);
    %sig0modelWu=sig0modelAll(5,:);
    sig0modelCM=sig0modelAll(8,:);
    
    clear sig0modelAll
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
        
    [surfFlag1,atmFrac]=makeSurfFlag(data,linInd);
    
    clear gasAttCloudMat
    %% Create field with reference sig0
    % RefFlag
    % 1 clear air
    % 2 interpolated
    % 3 model
    
    [refSig0,surfFlag,refFlag]=makeRefSig0(sig0measAtt,sig0modelCM,surfFlag1);
    
    % Find surfFlag values that have previously been clear air but are now
    % not
    surfFlag(surfFlag1~=0 & surfFlag==0 & any(data.FLAG==1,1))=1;
    
    %% 2 way path integrated attenuation from hydrometeors
    
    piaHydromet2=refSig0-sig0measAtt;
    piaHydromet2(surfFlag~=1)=nan;
    
    %% Separate warm and cold precip
    data.MELTING_LAYER(data.MELTING_LAYER<20)=10;
    data.MELTING_LAYER(data.MELTING_LAYER>19)=20;
    
    warmRefl=data.dbzMasked;
    warmRefl(data.MELTING_LAYER==20)=nan;
    coldRefl=data.dbzMasked;
    coldRefl(data.MELTING_LAYER==10)=nan;
    
    %% Make flag field that flags type of hydrometeor attenuation
    % 0 no attenuation
    % 1 liquid only
    % 2 mixed
    % 3 ice only
    
    warmFlag=any(~isnan(warmRefl),1);
    coldFlag=any(~isnan(coldRefl),1);
    
    attFlag=nan(size(data.time));
    attFlag(warmFlag & ~coldFlag)=1;
    attFlag(warmFlag & coldFlag)=2;
    attFlag(~warmFlag & coldFlag)=3;
    attFlag(surfFlag==0 | surfFlag==2)=0;
    
    %% Calculate two way ice attenuation
    
    % This equation comes from DOI: 10.1175/JTECH-D-18-0154.1 but I don't think
    % it applies to our very low reflectivities. Also, the units seem weird.
    %
    % coldReflLin=10.^(coldRefl./10);
    % iceSpecAtt=0.0325.*coldReflLin;
    %
    % iceAttAll=iceSpecAtt.*(data.range(2)-data.range(1))./1000;
    % piaIce2=sum(iceAttAll,1,'omitnan');
    
    clear coldRefl
    %% Calculate liquid attenuation
    
    piaLiq2=piaHydromet2;%-piaIce2;
    
    %% Calculate specific liquid attenuation with zPhi method
    
%     disp('Calculating specific liquid attenuation ...');
%     
%     specAttLiq=nan(size(data.DBZ));
%     
%     C1=4/(20*log10(exp(1)));
%     
%     b=nan(size(data.DBZ));
%     b(warmRefl<=-17)=b_drizz;
%     b(warmRefl>-17)=b_rain;
%     
%     meanB=mode(b,1);
%     
%     clear b
%     
%     cloudInds=find(attFlag==1);
%     
%     for ii=1:length(cloudInds)
%         dbzRay=warmRefl(:,cloudInds(ii));
%         cloudIndsRay=find(~isnan(dbzRay));
%         
%         if length(cloudIndsRay)>2
%             % Two way specific attenuation
%             % Z phi method
%             dbzLinB  = (10.^(0.1.*dbzRay)).^meanB(cloudInds(ii));
%             I0 = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay,cloudInds(ii))./1000,dbzLinB(cloudIndsRay));
%             CC = 10.^(0.1*meanB(cloudInds(ii))*piaLiq2(cloudInds(ii)))-1;
%             for mm = 1:length(cloudIndsRay)
%                 if mm < length(cloudIndsRay)
%                     Ir = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay(mm:end),cloudInds(ii))./1000,dbzLinB(cloudIndsRay(mm:end)));
%                 else
%                     Ir = 0;
%                 end
%                 specAttLiq(cloudIndsRay(mm),cloudInds(ii)) = (dbzLinB(cloudIndsRay(mm))*CC)/(I0+CC*Ir);
%             end
%         end
%     end
%     
%     specAttLiqNeg=find(specAttLiq<0);
%     specAttLiq(specAttLiqNeg)=nan;
%     
%     clear warmRefl cloudInds
%     %% Calculate LWC and RES with method taking Mie scattering into account
%     
%     disp('Calculating LWC, RES, and PID ...');
%     
%     PID=nan(size(data.DBZ));
%     LWC=nan(size(data.DBZ));
%     RES=nan(size(data.DBZ));
%     
%     piaInd=2*specAttLiq*(data.range(2)-data.range(1))./1000;
%     piaGate=cumsum(piaInd,1,'omitnan');
%     
%     clear piaInd
%     
%     piaGate(isnan(specAttLiq))=nan;
%       
%     % Attenuation corrected reflectivity
%     dbzCorr=data.dbzMasked+piaGate;

%% Plot

close all

sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);

    refSig0(upInds)=nan;
    data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY(upInds)=nan;
    
    sig0refMeas=nan(size(data.time));
    sig0refMeas(refFlag==1)=refSig0(refFlag==1);
    sig0refInt=nan(size(data.time));
    sig0refInt(refFlag==2)=refSig0(refFlag==2);
    sig0refMod=nan(size(data.time));
    sig0refMod(refFlag==3)=refSig0(refFlag==3);

timeInds=1:length(data.time);
   
    f1 = figure('Position',[200 500 1500 500],'DefaultAxesFontSize',12,'renderer','painters');
        
        hold on
        l0=plot(data.time(:,timeInds),sig0modelCM(:,timeInds),'-c','linewidth',2);
        l1=plot(data.time(:,timeInds),sig0measClear(:,timeInds),'-b','linewidth',1);
        l2=plot(data.time(:,timeInds),sig0measCloud(:,timeInds),'color',[0.5 0.5 0.5],'linewidth',0.5);
        l3=plot(data.time(:,timeInds),sig0refMeas(:,timeInds),'-r','linewidth',2);
        l4=plot(data.time(:,timeInds),sig0refInt(:,timeInds),'-','color',[0.5 0 1],'linewidth',2);
        l5=plot(data.time(:,timeInds),sig0refMod(:,timeInds),'-m','linewidth',2);
        ylabel('Sig0 (dB)');
        ylim([0 20]);
        
        yyaxis right
        l6=plot(data.time(:,timeInds),data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY(:,timeInds),'-k','linewidth',1);
        l7=plot(data.time(:,timeInds),piaLiq2(:,timeInds),'-g','linewidth',1);
        ylabel('Atten. (dB)');
        ylim([-5 15]);
        grid on
        set(gca,'YColor','k');
        
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        
        legend([l0 l1 l2 l3 l4 l5 l6 l7],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
            'sig0 ref meas','sig0 ref int','sig0 ref mod','2-way gas att','2-way PIA liq'},...
            'orientation','horizontal','location','north');
        title([datestr(data.time(1),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(end),'yyyy-mm-dd HH:MM:SS')])
                
        f1 = figure('Position',[200 500 1500 500],'DefaultAxesFontSize',12,'renderer','painters');
        
        colormap jet
        
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.DBZ(:,timeInds),'edgecolor','none');
        view(2);
        l1=plot(data.time(:,timeInds),data.altitude(:,timeInds)./1000,'-k','linewidth',2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([-0.5 ylimRefl]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        legend(l1,'Altitude');
        title('Reflectivity (dBZ)')
        
     %% Plot
%     
%     disp('Plotting convectivities ...');
%     
%     f1 = figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible',showPlot);
%     
%     s1=subplot(5,1,1);
%     
%     colormap jet
%     
%     hold on
%     surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([-35 25]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     colorbar
%     grid on
%     title('Reflectivity (dBZ)')
%     s1pos=s1.Position;
%     s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
%     
%     s2=subplot(5,1,2);
%     hold on
%     surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([0 1]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     colorbar
%     grid on
%     title('Reflectivity convectivity')
%     s2pos=s2.Position;
%     s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
% 
%     s3=subplot(5,1,3);
%     
%     hold on
%     surf(data.time,data.asl./1000,convBoth,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([0 0.7]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     colorbar
%     grid on
%     title('Reflectivity convectivity * velocity convectivity');
%     s3pos=s3.Position;
%     s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
% 
%     s4=subplot(5,1,4);
%     
%     hold on
%     surf(data.time,data.asl./1000,convVEL,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([0 1]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     colorbar
%     grid on
%     title('Velocity convectivity')
%     s4pos=s4.Position;
%     s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
% 
%     s5=subplot(5,1,5);
%     hold on
%     surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([-5 5]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     colorbar
%     grid on
%     title('Radial velocity')
%     s5pos=s5.Position;
%     s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
%     
%     linkaxes([s1 s2 s3 s4 s5],'xy');
%     
%     set(gcf,'PaperPositionMode','auto')
%     print(f1,[figdir,project,'_convectivity_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
%     
end