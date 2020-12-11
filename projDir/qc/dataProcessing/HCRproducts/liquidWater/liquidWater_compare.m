% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

b_drizz = 0.52; % Z<-17 dBZ
b_rain = 0.68; % Z>-17 dBZ
%alpha = 0.21;

ylimUpper=2;
adjustZeroMeter=350; % Assume melting layer is adjustZeroMeter below zero degree altitude

meltArea=2000; % melting layer +/- meltArea (meters) is considered in strat conv velocity algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sci/romatsch/liquidWaterHCR/'];
figdir='/home/romatsch/plots/HCR/liquidWater/compare/';

dataDir=['/run/media/romatsch/RSF0006/rsf/hcr/',project,'2hz/']; % HCR
indir=['/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/',project,'/']; % HSRL

% Loop through cases
casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/liquidWaterCompare_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Get HCR data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    % HCR
    data.DBZ = [];
    data.LDR=[];
    data.WIDTH=[];
    data.VEL_CORR=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    frq=ncread(fileList{1},'frequency');
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    data.freq=ncread(fileList{1},'frequency');
    
    %% Get HSRL data
    
    fileListHSRL=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if isempty(fileListHSRL)
        error('No data found.')
    end
    
    dataLid=[];
    
    %HSRL data
    dataLid.HSRL_Aerosol_Backscatter_Coefficient=[];
    dataLid.HSRL_Volume_Depolarization=[];
    dataLid.HSRL_Aerosol_Extinction_Coefficient=[];
    
    dataVarsLid=fieldnames(dataLid);
    
    % Load data
    dataLid=read_HCR(fileListHSRL,dataLid,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVarsLid)
        if ~isfield(dataLid,dataVarsLid{ii})
            dataVarsLid{ii}=[];
        end
    end
    
    dataLid.HSRL_Aerosol_Backscatter_Coefficient(dataLid.HSRL_Aerosol_Backscatter_Coefficient<0.5e-6)=nan;
    
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
    
    surfMask=nan(size(data.time));
    
    %sort out non nadir pointing
    surfMask(data.elevation>-85)=0;
    
    %sort out land
    surfMask(data.TOPO>0)=0;
    
    % sort out data from below 2500m altitude
    %surfMask(data.altitude<2500)=0;
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Calculate reflectivity sum inside and outside ocean surface to
    % distinguish clear air and cloud
    reflTemp=data.DBZ;
    
    % Remove bang
    reflTemp(data.FLAG==6)=nan;
    
    reflLin=10.^(reflTemp./10);
    reflOceanLin=nan(size(data.time));
    reflNoOceanLin=nan(size(data.time));
    
    for ii=1:length(data.time)
        if (~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-5)) & ~isnan(maxGate(ii))
            reflRay=reflLin(:,ii);
            reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
            reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
        end
    end
    
    % Remove data where reflectivity outside of ocean swath is more than
    % 0.8
    clearAir=find(reflNoOceanLin<=0.8);
    surfMask(clearAir)=2;
    surfMask(isnan(reflOceanLin))=0;
    
    % Find cloud data
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
    
    surfMask(find(any(~isnan(data.dbzMasked),1) & surfMask~=0 & surfMask~=2))=1;
    
    % Remove noise source cal, ant trans, and missing
    surfMask(find(any(data.FLAG>9,1) & surfMask~=0))=0;
    
    % Remove extinct
    surfMask(find(any(data.FLAG==3,1) & surfMask~=0))=0;
    
    if ~max(surfMask)==0
        %% Find melting layer and separate warm and cold precip
        
        findMelt=f_meltLayer(data,adjustZeroMeter);
        zeroInds=find(findMelt==0);
        oneInds=find(findMelt==1);
        twoInds=find(findMelt==2);
        threeInds=find(findMelt==3);
        
        findMelt(zeroInds)=nan;
        
        %meltType=sum(findMelt,1,'omitnan');
        
        meltInd=nan(size(data.time));
        warmRefl=nan(size(data.DBZ));
        coldRefl=nan(size(data.DBZ));
        
        for ii=1:length(meltInd)
            meltRay=findMelt(:,ii);
            goodMeltInds=find(~isnan(meltRay));
            if ~isempty(goodMeltInds)
                meltInd(ii)=max(goodMeltInds);
            else
                meltInd(ii)=1;
            end
            
            warmRefl(meltInd(ii):end,ii)=data.dbzMasked(meltInd(ii):end,ii);
            coldRefl(1:meltInd(ii)-1,ii)=data.dbzMasked(1:meltInd(ii)-1,ii);
        end
        
        %% Calculate two way ice attenuation
        coldReflLin=10.^(coldRefl./10);
        iceSpecAtt=0.0325.*coldReflLin;
        
        iceAttAll=iceSpecAtt.*(data.range(2)-data.range(1))./1000;
        iceAtt=sum(iceAttAll,1,'omitnan');
        
        %% One way and two way gaseous attenuation
        
        [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_atten(frq/1e+9,data);
        gasAttCloud2=2*gasAttCloud';
        
        %% Calculate clear air and cloudy ocean reflectivity
        
        reflSurfNoGas=data.DBZ(linInd);
        
        % Add gaseous attenuation back in
        reflSurf=reflSurfNoGas+gasAttCloud2;
        
        % Remove data with clouds
        dbzClear=nan(size(data.time));
        dbzCloud=nan(size(data.time));
        
        dbzClear(surfMask==2)=reflSurf(surfMask==2);
        dbzCloud(surfMask==1)=reflSurf(surfMask==1);
        
        % Add ice attenuation back in
        dbzCloudTot=dbzCloud+iceAtt;
        
        % Clear air ocean refl
        clearShort=dbzClear;
        clearShort(isnan(clearShort))=[];
        
        meanClearShort=movmedian(clearShort,100,'omitnan');
        meanClear=nan(size(data.time));
        meanClear(~isnan(dbzClear))=meanClearShort;
        
        meanClear(1)=meanClear(min(find(~isnan(meanClear))));
        
        for jj=2:length(meanClear)
            if isnan(meanClear(jj))
                meanClear(jj)=meanClear(jj-1);
            end
        end
        
        %% Calculate liquid attenuation
        
        attLiq=nan(size(data.time));
        specAtt=nan(size(data.DBZ));
        
        C1=4/(20*log10(exp(1)));
        
        b=nan(size(data.DBZ));
        b(warmRefl<=-17)=b_drizz;
        b(warmRefl>-17)=b_rain;
        
        meanB=median(b,1,'omitnan');
        
        cloudInds=find(surfMask==1);
        
        for ii=1:length(cloudInds)
            dbzRay=warmRefl(:,cloudInds(ii));
            cloudIndsRay=find(~isnan(dbzRay));
            
            if length(cloudIndsRay)>2
                
                attDiff=meanClear(cloudInds(ii))-dbzCloudTot(cloudInds(ii));
                if attDiff>=0
                    attLiq(cloudInds(ii))=attDiff;
                else
                    attLiq(cloudInds(ii))=0;
                end
                
                % Two way specific attenuation
                % Z phi method
                dbzLinB  = (10.^(0.1.*dbzRay)).^meanB(cloudInds(ii));
                I0 = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay,cloudInds(ii))./1000,dbzLinB(cloudIndsRay));
                CC = 10.^(0.1*meanB(cloudInds(ii))*attLiq(cloudInds(ii)))-1;
                for mm = 1:length(cloudIndsRay)
                    if mm < length(cloudIndsRay)
                        Ir = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay(mm:end),cloudInds(ii))./1000,dbzLinB(cloudIndsRay(mm:end)));
                    else
                        Ir = 0;
                    end
                    specAtt(cloudIndsRay(mm),cloudInds(ii)) = (dbzLinB(cloudIndsRay(mm))*CC)/(I0+CC*Ir);
                end
            end
        end
        
        alpha=1./(4.792-3.63e-2*data.TEMP-1.897e-4*data.TEMP.^2);
        LWC=specAtt.*alpha;
        
        %% From HSRL
        
        % Calculate attenuation correction
        
        Z_95_lin=10.^(data.DBZ*0.1);
        Z_95_lin(data.DBZ < -200)=0.;
        
        wt_coef=zeros(size(data.DBZ));
        wt_exp=zeros(size(data.DBZ));
        wt_coef(data.DBZ < - 20)=20.;
        wt_exp(data.DBZ < - 20)=0.52;
        wt_coef(-20 <data.DBZ <-15 )=1.73;
        wt_exp(-20 <data.DBZ < -15 )=0.15;
        wt_coef(data.DBZ > -15)=0.22;
        wt_exp(data.DBZ > -15)=0.68;
        
        att_cumul=2.*0.0192*cumsum((wt_coef.*Z_95_lin.^wt_exp),1,'omitnan');
        att_cumul(data.DBZ < -200)=NaN;
        
        dBZ_cor=data.DBZ+att_cumul;
        Z_95_lin_cor=10.^(dBZ_cor*0.1);
        
        %% Censor data
        
        cen_ind=dataLid.HSRL_Aerosol_Backscatter_Coefficient  <= 0;
        
        dBZ = data.DBZ;
        dBZ(cen_ind)= NaN;
        beta = dataLid.HSRL_Aerosol_Backscatter_Coefficient;
        beta(cen_ind) = NaN;
        
        %% Calculate RES and LWC
        
        radar_lidar= dBZ_cor-10.*log10(beta);
        radar_lidar_lin= 10.^(radar_lidar*0.1);
        
        retrieved_RES=8.3478*(radar_lidar_lin).^0.3102; % in microns
        ref_normRES_retrieved=Z_95_lin_cor./((1e-3*retrieved_RES).^3); % retrieved Cld.RES is in microns and it is converted to mm by the fctor 1e-3
        dsdretrieved_LWC=5.238e-4*ref_normRES_retrieved+2.212e-4;
        
        dsdretrieved_LWC(data.DBZ>0)=nan;
        
        LWCdiff=LWC-dsdretrieved_LWC;
        LWCdiff(LWC==0)=nan;
        
        %% Plot liquid attenuation
        close all
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
        
        s1=subplot(3,1,1);
        hold on
        l1=plot(data.time,dbzClear,'-b','linewidth',1);
        l2=plot(data.time,dbzCloud,'color',[0.5 0.5 0.5],'linewidth',0.5);
        l3=plot(data.time,meanClear,'-r','linewidth',2);
        ylabel('Refl. (dBZ)');
        ylim([40 60]);
        
        yyaxis right
        l4=plot(data.time,gasAttCloud2,'-k','linewidth',1);
        l5=plot(data.time,attLiq,'-g','linewidth',1);
        l6=plot(data.time,iceAtt*10,'-m','linewidth',1);
        ylabel('Atten. (dB)');
        ylim([-5 15]);
        grid on
        set(gca,'YColor','k');
        
        xlim([data.time(1),data.time(end)]);
        
        legend([l1 l3 l4 l5 l6],{'Refl. measured','Refl. used','2-way gaseous atten.','2-way liquid atten.','2-way ice atten. * 10'},...
            'orientation','horizontal','location','south');
        title([datestr(data.time(1)),' to ',datestr(data.time(end))])
        s1pos=s1.Position;
        
        s2=subplot(3,1,2);
        
        colormap jet
        
        hold on
        surf(data.time,data.asl./1000,data.dbzMasked,'edgecolor','none');
        view(2);
        scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'k','filled');
        scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'c','filled');
        scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'b','filled');
        scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
        ax = gca;
        ax.SortMethod = 'childorder';
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        colorbar
        grid on
        title('Reflectivity (dBZ)')
        s2pos=s2.Position;
        s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
        
        s3=subplot(3,1,3);
        
        colmap=jet;
        colmap=cat(1,[1 0 1],colmap);
        
        hold on
        surf(data.time,data.asl./1000,LWC,'edgecolor','none');
        view(2);
        colormap(s3,colmap)
        ylabel('Altitude (km)');
        caxis([0 2]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        colorbar
        grid on
        title('Liquid water content (g m^{-3})')
        s3pos=s3.Position;
        s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_lwcHCR_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        %% Plot liquid attenuation from HCR, HSRL, and difference
        
        f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
        
        s1=subplot(3,1,1);
        
        colmap=jet;
        colmap=cat(1,[1 0 1],colmap);
        
        hold on
        surf(data.time,data.asl./1000,LWC,'edgecolor','none');
        view(2);
        colormap(s1,colmap)
        ylabel('Altitude (km)');
        caxis([0 2]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        colorbar
        grid on
        title('Liquid water content HCR (g m^{-3})')
        
        s2=subplot(3,1,2);
        hold on
        surf(data.time,data.asl./1000,dsdretrieved_LWC,'edgecolor','none');
        view(2);
        colormap(s2,colmap)
        ylabel('Altitude (km)');
        caxis([0 2]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        colorbar
        grid on
        title('Liquid water content HCRandHSRL (g m^{-3})')
        
        s3=subplot(3,1,3);
        
        colmap=jet;
        colmap=cat(1,[1 0 1],colmap);
        
        hold on
        surf(data.time,data.asl./1000,LWCdiff,'edgecolor','none');
        view(2);
        colormap(s3,jet)
        ylabel('Altitude (km)');
        caxis([-2 2]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        colorbar
        grid on
        title('Liquid water content HCR-HCRandHSRL (g m^{-3})')
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_lwcCompare_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    end
end
