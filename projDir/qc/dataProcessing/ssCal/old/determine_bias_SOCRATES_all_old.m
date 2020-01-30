% Plots flight and HCR info to be used in determination of surface
% reflectivity.  Determines bias according to Rilling technique.
% Uses a calibration input file, with attenuations pre-computed
% by Romatschke.
%
% Creates multi-panel plots over a elected, contiguous, number of
% data files.  Panels include:
%
%   -- scan rotation angle (rotation is in a plane perpendicular
%                           to flight direction)
%   -- max reflectivity (when pointing down, only)
%   -- range to max reflectivity (when pointing down)
%   -- Altitude
%   -- Pitch
%
% all parameters are plotted vs time, where time is in seconds from
% the start of the first file, or seconds from the start of the day.
%
% Note that for sea surface scans, we operate in vertical receive, only.

clear all;
close all;

makeFigs=0;

addpath('/h/eol/romatsch/git/private/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/git/private/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/git/private/process_HCR/oceanScans/subCodes/');

figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/socrates/';

indir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/';
sondedir='/h/eol/romatsch/data/hcrCalib/soundings/socrates_drops/';
modeldir='/scr/sci/romatsch/data/reanalysis/ecmwf/socrates/soundIn/';

filedir='/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/';
infile='cal_SOCRATES_slow.txt';

frq = 9.4406e10;
%% Slow scans

sfcWdspdSlow=[];
sfcWddirSlow=[];
lbBiasSlow=[];
ituBiasSlow=[];
hdgSlow=[];
altSlow=[];
attLbSlow=[];
attITUSlow=[];
bias_lbITU_vec_slow=[];

% read file with cases
caseList = readtable([filedir,infile],'Delimiter', ' ');
%convert to cell so each case has one cell
casesIn=table2array(caseList(:,1));
numCases=unique(casesIn);
uniqueCasesSlow=cell(size(numCases,1),1);

PLTslow=cell(size(numCases,1),1);

for ii=1:size(numCases,1)
    caseInd=find(casesIn==ii);
    uniqueCasesSlow{ii}=caseList(caseInd,:);
end


for ii=1:size(uniqueCasesSlow,1); %Loop through all cases
    disp(['Case ',num2str(ii),' of ',num2str(size(uniqueCasesSlow,1))]);
       
    if strcmp(uniqueCasesSlow{ii,1}.soundg{1},'nan')
        uniqueCasesSlow{ii,1}.soundg{1}=str2num(uniqueCasesSlow{ii,1}.soundg{1});
    end
    
    if ~isnan(uniqueCasesSlow{ii,1}.soundg{1})
        fileType=uniqueCasesSlow{ii,1}.soundg{1}(end-2:end);
        if strcmp(fileType,'eol')
            soundFile=1;
        elseif strcmp(fileType,'.nc')
            soundFile=2;
        else
            disp('Wrong sounding or model file extension.');
        end
        
        if soundFile==1
            % calculate attenuation from sounding
            [attLiebe,attITU,sfcWdspd,sfcWddir]=f_atten_layers_sfcWnd([sondedir,uniqueCasesSlow{ii,1}.soundg{1}],frq/1e+9);
        elseif soundFile==2
            % calculate attenuation from ecmwf model data
            %get lat lon alt
            oneCaseFileList=importdata([filedir,uniqueCasesSlow{ii,1}.casefiles{1}]);
            
            indata=[indir,oneCaseFileList{1}];
            
            lon=ncread(indata,'longitude');
            lat=ncread(indata,'latitude');
            alt=ncread(indata,'altitude');
            
            %model
            [attLiebe,attITU,sfcWdspd,sfcWddir]=f_atten_layers_sfcWnd_ecmwf([modeldir,uniqueCasesSlow{ii,1}.soundg{1}],...
                frq/1e+9,nanmean(lon),nanmean(lat),nanmean(alt));
        end
    else
        attLiebe=nan;
        attITU=nan;
        sfcWdspd=nan;
        sfcWddir=nan;
    end
    
    sfcWdspdSlow=[sfcWdspdSlow sfcWdspd];
    sfcWddirSlow=[sfcWddirSlow sfcWdspd];
        
     
    [PLT lbBias ituBias avg_hdg avg_alt bias_lbITU_vec_slow]=f_determine_bias(uniqueCasesSlow{ii,1},filedir,indir,...
        attITU,attLiebe,frq,bias_lbITU_vec_slow);
    
    PLTslow{ii}=PLT;
    lbBiasSlow=[lbBiasSlow lbBias];
    ituBiasSlow=[ituBiasSlow ituBias];
    hdgSlow=[hdgSlow avg_hdg];
    altSlow=[altSlow avg_alt];
    attLbSlow=[attLbSlow attLiebe];
    attITUSlow=[attITUSlow attITU];
end

%% Fast scans

sfcWdspdFast=[];
sfcWddirFast=[];
lbBiasFast=[];
ituBiasFast=[];
hdgFast=[];
altFast=[];
attLbFast=[];
attITUFast=[];
bias_lbITU_vec_fast=[];

infile='cal_SOCRATES_fast.txt';

% read file with cases
caseList = readtable([filedir,infile],'Delimiter','space');
%convert to cell so each case has one cell
casesIn=table2array(caseList(:,1));
numCases=unique(casesIn);
uniqueCasesFast=cell(size(numCases,1),1);

PLTfast=cell(size(numCases,1),1);

for ii=1:size(numCases,1)
    caseInd=find(casesIn==ii);
    uniqueCasesFast{ii}=caseList(caseInd,:);
end


for ii=1:size(uniqueCasesFast,1); %Loop through all cases
    disp(['Case ',num2str(ii),' of ',num2str(size(uniqueCasesFast,1))]);
    
    if strcmp(uniqueCasesFast{ii,1}.soundg{1},'nan')
        uniqueCasesFast{ii,1}.soundg{1}=str2num(uniqueCasesFast{ii,1}.soundg{1});
    end
    
    if ~isnan(uniqueCasesFast{ii,1}.soundg{1})
        fileType=uniqueCasesFast{ii,1}.soundg{1}(end-2:end);
        if strcmp(fileType,'eol')
            soundFile=1;
        elseif strcmp(fileType,'.nc')
            soundFile=2;
        else
            disp('Wrong sounding or model file extension.');
        end
        
        if soundFile==1
            % calculate attenuation from sounding
            [attLiebe,attITU,sfcWdspd,sfcWddir]=f_atten_layers_sfcWnd([sondedir,uniqueCasesFast{ii,1}.soundg{1}],frq/1e+9);
        elseif soundFile==2
            % calculate attenuation from ecmwf model data
            %get lat lon alt
            oneCaseFileList=importdata([filedir,uniqueCasesFast{ii,1}.casefiles{1}]);
            
            indata=[indir,oneCaseFileList{1}];
            
            lon=ncread(indata,'longitude');
            lat=ncread(indata,'latitude');
            alt=ncread(indata,'altitude');
            
            %model
            [attLiebe,attITU,sfcWdspd,sfcWddir]=f_atten_layers_sfcWnd_ecmwf([modeldir,uniqueCasesFast{ii,1}.soundg{1}],...
                frq/1e+9,nanmean(lon),nanmean(lat),nanmean(alt));
        end
    else
        attLiebe=nan;
        attITU=nan;
        sfcWdspd=nan;
        sfcWddir=nan;
    end
    
    sfcWdspdFast=[sfcWdspdFast sfcWdspd];
    sfcWddirFast=[sfcWddirFast sfcWdspd];
    
    [PLT lbBias ituBias avg_hdg avg_alt bias_lbITU_vec_fast]=f_determine_bias(uniqueCasesFast{ii,1},filedir,indir,...
        attITU,attLiebe,frq,bias_lbITU_vec_fast);
    
    PLTfast{ii}=PLT;
    lbBiasFast=[lbBiasFast lbBias];
    ituBiasFast=[ituBiasFast ituBias];
    hdgFast=[hdgFast avg_hdg];
    altFast=[altFast avg_alt];
    attLbFast=[attLbFast attLiebe];
    attITUFast=[attITUFast attITU];
end

%% Create table

%fast
fastTable=cell2table(cell(size(uniqueCasesFast,1),10), ...
    'VariableNames', {'StartDate','EndDate','OneWayAttLiebe','OneWayAttITU','BiasLiebe','BiasITU', ...
    'MeanHdg','Sounding','SfcWindSpd','SfcWindDir'});
for ii=1:size(fastTable,1)
    fastTable(ii,1)=uniqueCasesFast{ii,1}.timest(1);
    fastTable(ii,2)=uniqueCasesFast{ii,1}.timend(end);
    soundingTemp=uniqueCasesFast{ii,1}.soundg(1);
    if isnan(soundingTemp{:})
        fastTable(ii,8)={'nan'};
    else
        fastTable(ii,8)=soundingTemp;
    end
end

fastTable.OneWayAttLiebe=attLbFast';
fastTable.OneWayAttITU=attITUFast';
fastTable.BiasLiebe=lbBiasFast';
fastTable.BiasITU=ituBiasFast';
fastTable.MeanHdg=hdgFast';
fastTable.SfcWindSpd=sfcWdspdFast';
fastTable.SfcWindDir=sfcWddirFast';

%slow
slowTable=cell2table(cell(size(uniqueCasesSlow,1),10), ...
    'VariableNames', {'StartDate','EndDate','OneWayAttLiebe','OneWayAttITU','BiasLiebe','BiasITU', ...
    'MeanHdg','Sounding','SfcWindSpd','SfcWindDir'});
for ii=1:size(slowTable,1)
    slowTable(ii,1)=uniqueCasesSlow{ii,1}.timest(1);
    slowTable(ii,2)=uniqueCasesSlow{ii,1}.timend(end);
    soundingTemp=uniqueCasesSlow{ii,1}.soundg(1);
    if isnan(soundingTemp{:})
        slowTable(ii,8)={'nan'};
    else
        slowTable(ii,8)=soundingTemp;
    end
end

slowTable.OneWayAttLiebe=attLbSlow';
slowTable.OneWayAttITU=attITUSlow';
slowTable.BiasLiebe=lbBiasSlow';
slowTable.BiasITU=ituBiasSlow';
slowTable.MeanHdg=hdgSlow';
slowTable.SfcWindSpd=sfcWdspdSlow';
slowTable.SfcWindDir=sfcWddirSlow';

%% Plot time series
if makeFigs
    
    close all
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        f_plot_series(PLT,[figdir 'timeSeries/SOCRATES_fast_series_' uniqueCasesFast{ii}.timest{1}]);
    end
    close all
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        f_plot_series(PLT,[figdir 'timeSeries/SOCRATES_slow_series_' uniqueCasesSlow{ii}.timest{1}]);
    end
    
    %% Plot sig0 for individual cases
    close all
    xlims=[0 25];
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        titleIn=['Bias-corrected Sigma0 vs Incidence Angle ' uniqueCasesFast{ii}.timest{1}];
        if min(isnan(PLT.sig0))==0
            f_plot_sig0(PLT,[figdir 'sigma0/SOCRATES_fast_sigma0_' uniqueCasesFast{ii}.timest{1}],...
                xlims,titleIn,sfcWdspdFast(ii),sfcWdspdFast(ii),hdgFast(ii),altFast(ii),lbBiasFast(ii),attLbFast(ii));
        end
    end
    
    close all
    xlims=[5 15];
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        titleIn=['Bias-corrected Sigma0 vs Incidence Angle ' uniqueCasesSlow{ii}.timest{1}];
        if min(isnan(PLT.sig0))==0
            f_plot_sig0(PLT,[figdir 'sigma0/SOCRATES_slow_sigma0_' uniqueCasesSlow{ii}.timest{1}],...
                xlims,titleIn,sfcWdspdSlow(ii),sfcWdspdSlow(ii),hdgSlow(ii),altSlow(ii),lbBiasSlow(ii),attLbSlow(ii));
        end
    end
      %% Plot sig0noBias for individual cases
    close all
    xlims=[0 25];
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        titleIn=['Non-corrected Sigma0 vs Incidence Angle ' uniqueCasesFast{ii}.timest{1}];
        if min(isnan(PLT.sig0))==0
            f_plot_sig0noBias(PLT,[figdir 'sigma0noBias/SOCRATES_fast_sigma0noBias_' uniqueCasesFast{ii}.timest{1}],...
                xlims,titleIn);
        end
    end
    
    close all
    xlims=[5 15];
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        titleIn=['Non-corrected Sigma0 vs Incidence Angle ' uniqueCasesFast{ii}.timest{1}];
        if min(isnan(PLT.sig0))==0
            f_plot_sig0noBias(PLT,[figdir 'sigma0noBias/SOCRATES_slow_sigma0noBias_' uniqueCasesFast{ii}.timest{1}],...
                xlims,titleIn);
        end
    end
    %% Plot reflectivity for individual cases
    close all
    xlims=[0 25];
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        titleIn=['Reflectivity vs Incidence Angle ' uniqueCasesFast{ii}.timest{1}];
        if min(isnan(PLT.refl))==0
            f_plot_refl(PLT,[figdir 'reflectivity/SOCRATES_fast_refl_' uniqueCasesFast{ii}.timest{1}],xlims,titleIn);
        end
    end
    
    close all
    xlims=[5 15];
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        titleIn=['Reflectivity vs Incidence Angle ' uniqueCasesSlow{ii}.timest{1}];
        if min(isnan(PLT.refl))==0
            f_plot_refl(PLT,[figdir 'reflectivity/SOCRATES_slow_refl_' uniqueCasesSlow{ii}.timest{1}],xlims,titleIn);
        end
    end
    
end
    %% Plot all sig0 vs elev angle and sig0vs wind speed
    close all
    
    f1=figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    f3 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTfast,1)
        if min(isnan(PLTfast{ii,1}.sig0))==0
            
            wdspFast=nan(size(PLTfast{ii,1}.sig0));
            wdspFast(:,:)=sfcWdspdFast(ii);
            
            sig0smallerFast=PLTfast{ii,1}.sig0(find(PLTfast{ii,1}.rota<=180));
            sig0noBiassmallerFast=PLTfast{ii,1}.sig0noBias(find(PLTfast{ii,1}.rota<=180));
            elevSmallerFast=PLTfast{ii,1}.elev(find(PLTfast{ii,1}.rota<=180));
            wdSmallerFast=wdspFast(find(PLTfast{ii,1}.rota<=180));
            
            sig0biggerFast=PLTfast{ii,1}.sig0(find(PLTfast{ii,1}.rota>180));
            sig0noBiasbiggerFast=PLTfast{ii,1}.sig0noBias(find(PLTfast{ii,1}.rota>180));
            elevBiggerFast=PLTfast{ii,1}.elev(find(PLTfast{ii,1}.rota>180));
            wdBiggerFast=wdspFast(find(PLTfast{ii,1}.rota>180));
            
            figure(f1)
            plot(elevSmallerFast,sig0smallerFast,'r');
            plot(elevBiggerFast,sig0biggerFast,'b');
            
            figure(f3)
            plot(elevSmallerFast,sig0noBiassmallerFast,'r');
            plot(elevBiggerFast,sig0noBiasbiggerFast,'b');
            
            elevIndsSF=find(elevSmallerFast>9.7 & elevSmallerFast<10.3);
            elevIndsBF=find(elevBiggerFast>9.7 & elevBiggerFast<10.3);
            
            figure(f2)
            plot(wdSmallerFast(elevIndsSF),sig0smallerFast(elevIndsSF),'or');
            plot(wdBiggerFast(elevIndsBF),sig0biggerFast(elevIndsBF),'ob');
        end
    end
    for ii=1:size(PLTslow,1)
        if min(isnan(PLTslow{ii,1}.sig0))==0
            
            wdspSlow=nan(size(PLTslow{ii,1}.sig0));
            wdspSlow(:,:)=sfcWdspdSlow(ii);
            
            sig0smallerSlow=PLTslow{ii,1}.sig0(find(PLTslow{ii,1}.rota<=180));
            sig0noBiassmallerSlow=PLTslow{ii,1}.sig0noBias(find(PLTslow{ii,1}.rota<=180));
            elevSmallerSlow=PLTslow{ii,1}.elev(find(PLTslow{ii,1}.rota<=180));
            wdSmallerSlow=wdspFast(find(PLTslow{ii,1}.rota<=180));
            
            sig0biggerSlow=PLTslow{ii,1}.sig0(find(PLTslow{ii,1}.rota>180));
            sig0noBiasbiggerSlow=PLTslow{ii,1}.sig0noBias(find(PLTslow{ii,1}.rota>180));
            elevBiggerSlow=PLTslow{ii,1}.elev(find(PLTslow{ii,1}.rota>180));
            wdBiggerSlow=wdspSlow(find(PLTslow{ii,1}.rota>180));
            
            figure(f1)
            plot(elevSmallerSlow,sig0smallerSlow,'m');
            plot(elevBiggerSlow,sig0biggerSlow,'c');
            
            figure(f3)
            plot(elevSmallerSlow,sig0noBiassmallerSlow,'m');
            plot(elevBiggerSlow,sig0noBiasbiggerSlow,'c');
            
            elevIndsSS=find(elevSmallerSlow>9.7 & elevSmallerSlow<10.3);
            elevIndsBS=find(elevBiggerSlow>9.7 & elevBiggerSlow<10.3);
            
            figure(f2)
            plot(wdSmallerSlow(elevIndsSS),sig0smallerSlow(elevIndsSS),'om');
            plot(wdBiggerSlow(elevIndsBS),sig0biggerSlow(elevIndsBS),'oc');
        end
    end
    
    figure(f1)
    set(gca,'XLim',[0 25]);
    set(gca,'YLim',[-20,20]);
    xlabel('Incidence Angle, degrees');
    ylabel('Normalized Radar Cross Section, dB');
    title(['Bias-corrected Sigma0 vs Incidence Angle']);
    legend('Rot angle <= 180 deg fast','Rot angle > 180 deg fast','Rot angle <= 180 deg slow','Rot angle > 180 deg slow','location','southwest');
    plabel = [figdir 'SOCRATES_sigma0_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    figure(f3)
    set(gca,'XLim',[0 25]);
    set(gca,'YLim',[-20,20]);
    xlabel('Incidence Angle, degrees');
    ylabel('Raw normalized Radar Cross Section, dB');
    title(['Non-corrected Sigma0 vs Incidence Angle']);
    legend('Rot angle <= 180 deg fast','Rot angle > 180 deg fast','location','southwest');
    plabel = [figdir 'SOCRATES_sigma0noBias_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    figure(f2)
    set(gca,'YLim',[-5,10]);
    xlabel('Surface wind speed, m/s');
    ylabel('Normalized Radar Cross Section, dB');
    title(['Bias-corrected Sigma0 vs Surface wind speed']);
    legend('Rot angle <= 180 deg fast','Rot angle > 180 deg fast','Rot angle <= 180 deg slow','Rot angle > 180 deg slow','location','southwest');
    plabel = [figdir 'SOCRATES_sigma0_windspd_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias vs wind speed
    close all
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(sfcWdspdFast,lbBiasFast,'o', 'MarkerFaceColor', 'b');
    plot(sfcWdspdSlow,lbBiasSlow, 'oc', 'MarkerFaceColor', 'c');
    set(gca,'XLim',[0 20]);
    set(gca,'YLim',[0,4]);
    xlabel('Surface wind speed, m/s');
    ylabel('LiebeBias, dB');
    title(['Liebe bias vs Surface wind speed']);
    legend('Fast','Slow','location','southeast');
    plabel = [figdir 'SOCRATES_bias_windspd_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias vs wind angle
    
    hdgFast(hdgFast<0)=hdgFast(hdgFast<0)+360;
    hdgSlow(hdgSlow<0)=hdgSlow(hdgSlow<0)+360;
    
    normDegF = mod(hdgFast-sfcWddirFast,360);
    hdgMinusDirF=min(360-normDegF, normDegF);
    normDegS = mod(hdgSlow-sfcWddirSlow,360);
    hdgMinusDirS=min(360-normDegS, normDegS);
    
    close all
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(hdgMinusDirF,lbBiasFast,'o', 'MarkerFaceColor', 'b');
    plot(hdgMinusDirS,lbBiasSlow,'oc', 'MarkerFaceColor', 'c');
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0,4]);
    xlabel('abs(heading - surface wind dir), degree');
    ylabel('LiebeBias, dB');
    title(['Liebe bias vs heading minus wind direction']);
    legend('Fast','Slow','location','southwest');
    plabel = [figdir 'SOCRATES_bias_winddir_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% To smooth, put data in bins
    
    edgesFast=(0:0.5:25);
    
    smElevF=[];
    smReflF=[];
    smSig0F=[];
    smRotaF=[];
    smSig0noBiasF=[];
    
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        if isnan(PLT.sig0)
            PLT.sig0=nan(size(PLT.elev));
            PLT.sig0noBias=nan(size(PLT.elev));
        end
        elRefSigRot=cat(2,PLT.elev,PLT.refl,PLT.sig0,PLT.rota,PLT.sig0noBias);
        binInds=discretize(PLT.elev,edgesFast);
        
        binData=nan(length(edgesFast-1),5);
        
        for jj=1:length(edgesFast-1)
            rowInds=find(binInds==jj);
            data=elRefSigRot(rowInds,:);
            binData(jj,1)=mean(data(:,1));
            [binData(jj,2),~,~]=dB_meanStd(data(:,2));
            [binData(jj,3),~,~]=dB_meanStd(data(:,3));
            binData(jj,4)=mean(data(:,4));
            [binData(jj,5),~,~]=dB_meanStd(data(:,5));
        end
        
        smElevF=cat(2,smElevF,binData(:,1));
        smReflF=cat(2,smReflF,binData(:,2));
        smSig0F=cat(2,smSig0F,binData(:,3));
        smRotaF=cat(2,smRotaF,binData(:,4));
        smSig0noBiasF=cat(2,smSig0noBiasF,binData(:,5));
    end
    
    edgesSlow=(0:0.5:25);
    
    smElevS=[];
    smReflS=[];
    smSig0S=[];
    smRotaS=[];
    smSig0noBiasS=[];
    
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        if isnan(PLT.sig0)
            PLT.sig0=nan(size(PLT.elev));
            PLT.sig0noBias=nan(size(PLT.elev));
        end
        elRefSigRot=cat(2,PLT.elev,PLT.refl,PLT.sig0,PLT.rota,PLT.sig0noBias);
        binInds=discretize(PLT.elev,edgesSlow);
        
        binData=nan(length(edgesSlow-1),5);
        
        for jj=1:length(edgesSlow-1)
            rowInds=find(binInds==jj);
            data=elRefSigRot(rowInds,:);
            binData(jj,1)=mean(data(:,1));
            [binData(jj,2),~,~]=dB_meanStd(data(:,2));
            [binData(jj,3),~,~]=dB_meanStd(data(:,3));
            binData(jj,4)=mean(data(:,4));
            [binData(jj,5),~,~]=dB_meanStd(data(:,5));
        end
        
        smElevS=cat(2,smElevS,binData(:,1));
        smReflS=cat(2,smReflS,binData(:,2));
        smSig0S=cat(2,smSig0S,binData(:,3));
        smRotaS=cat(2,smRotaS,binData(:,4));
        smSig0noBiasS=cat(2,smSig0noBiasS,binData(:,5));
    end
    
    %% Plot smoothed data reflectivity
    close all
    
    %fast scans
    legendF=[];
    colmapFast;
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTfast,1)
        legendF{end+1}=uniqueCasesFast{ii}.timest{1};
        plot(smElevF(:,ii),smReflF(:,ii),'Color',colMapF(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[0 35]);
    set(gca,'YLim',[25,55]);
    xlabel('Incidence Angle, degrees');
    ylabel('Reflectivity, dB');
    title(['Reflectivity vs Incidence Angle smoothed']);
    legend(legendF,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_fast_refl_all_smooth'],'-dpng','-r0')
    
    % slow scans
    legendS=[];
    colmapSlow;
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTslow,1)
        legendS{end+1}=uniqueCasesSlow{ii}.timest{1};
        plot(smElevS(:,ii),smReflS(:,ii),'Color',colMapS(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[5 20]);
    set(gca,'YLim',[25,55]);
    xlabel('Incidence Angle, degrees');
    ylabel('Reflectivity, dB');
    title(['Reflectivity vs Incidence Angle smoothed']);
    legend(legendS,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_slow_refl_all_smooth'],'-dpng','-r0')
    
    %% Plot smoothed data sig0
    close all
    
    %fast scans
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTfast,1)
        plot(smElevF(:,ii),smSig0F(:,ii),'Color',colMapF(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[0 35]);
    set(gca,'YLim',[-10,15]);
    xlabel('Incidence Angle, degrees');
    ylabel('Normalized Radar Cross Section, dB');
    title(['Bias-corrected Sigma0 vs Incidence Angle, smoothed']);
    legend(legendF,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_fast_sigma0_all_smooth'],'-dpng','-r0')
    
    % slow scans
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTslow,1)
        plot(smElevS(:,ii),smSig0S(:,ii),'Color',colMapS(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[5 20]);
    set(gca,'YLim',[-10,15]);
    xlabel('Incidence Angle, degrees');
    ylabel('Normalized Radar Cross Section, dB');
    title(['Bias-corrected Sigma0 vs Incidence Angle, smoothed']);
    legend(legendS,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_slow_sigma0_all_smooth'],'-dpng','-r0')
    
    %% Plot smoothed data sig0noBias
    close all
    
    %fast scans
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTfast,1)
        plot(smElevF(:,ii),smSig0noBiasF(:,ii),'Color',colMapF(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[0 35]);
    set(gca,'YLim',[-10,15]);
    xlabel('Incidence Angle, degrees');
    ylabel('Raw normalized Radar Cross Section, dB');
    title(['Non-corrected Sigma0 vs Incidence Angle, smoothed']);
    legend(legendF,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_fast_sigma0noBias_all_smooth'],'-dpng','-r0')
    
    % slow scans
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    
    for ii=1:size(PLTslow,1)
        plot(smElevS(:,ii),smSig0noBiasS(:,ii),'Color',colMapS(ii,:),'linewidth',1.5);
    end
    
    set(gca,'XLim',[5 20]);
    set(gca,'YLim',[-10,15]);
    xlabel('Incidence Angle, degrees');
    ylabel('Raw normalized Radar Cross Section, dB');
    title(['Non-corrected Sigma0 vs Incidence Angle, smoothed']);
    legend(legendS,'location','northeast','Interpreter','none');
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_slow_sigma0noBias_all_smooth'],'-dpng','-r0')
    %% 
    %Calculate mean and std over all data, not means of means
    lbITUmeanS=nan(1,2);
    lbITUstdUS=nan(1,2);
    lbITUstdDS=nan(1,2);
    [lbITUmeanS(1),lbITUstdUS(1),lbITUstdDS(1)]=dB_meanStd(bias_lbITU_vec_slow(:,1));
    [lbITUmeanS(2),lbITUstdUS(2),lbITUstdDS(2)]=dB_meanStd(bias_lbITU_vec_slow(:,2));
   
    lbITUmeanF=nan(1,2);
    lbITUstdUF=nan(1,2);
    lbITUstdDF=nan(1,2);
    [lbITUmeanF(1),lbITUstdUF(1),lbITUstdDF(1)]=dB_meanStd(bias_lbITU_vec_fast(:,1));
    [lbITUmeanF(2),lbITUstdUF(2),lbITUstdDF(2)]=dB_meanStd(bias_lbITU_vec_fast(:,2));
    
    lbITUmeanAll=nan(1,2);
    lbITUstdUAll=nan(1,2);
    lbITUstdDAll=nan(1,2);
    [lbITUmeanAll(1),lbITUstdUAll(1),lbITUstdDAll(1)]=dB_meanStd(cat(1,bias_lbITU_vec_slow(:,1),bias_lbITU_vec_fast(:,1)));
    [lbITUmeanAll(2),lbITUstdUAll(2),lbITUstdDAll(2)]=dB_meanStd(cat(1,bias_lbITU_vec_slow(:,2),bias_lbITU_vec_fast(:,2)));
    
    disp(['Fast Liebe M: ',num2str(lbITUmeanF(1)),' + ',num2str(lbITUstdUF(1)),' - ',num2str(lbITUstdDF(1))]);
    disp(['Fast ITU M: ',num2str(lbITUmeanF(2)),' + ',num2str(lbITUstdUF(2)),' - ',num2str(lbITUstdDF(2))]);
    
    disp(['Slow Liebe M: ',num2str(lbITUmeanS(1)),' + ',num2str(lbITUstdUS(1)),' - ',num2str(lbITUstdDS(1))]);
    disp(['Slow ITU M: ',num2str(lbITUmeanS(2)),' + ',num2str(lbITUstdUS(2)),' - ',num2str(lbITUstdDS(2))]);
        
    disp(['All Liebe M: ',num2str(lbITUmeanAll(1)),' + ',num2str(lbITUstdUAll(1)),' - ',num2str(lbITUstdDAll(1))]);
    disp(['All ITU M: ',num2str(lbITUmeanAll(2)),' + ',num2str(lbITUstdUAll(2)),' - ',num2str(lbITUstdDAll(2))]);
    
    [biasLmeanF,biasLstdUF,biasLstdDF]=dB_meanStd(fastTable.BiasLiebe);
    [biasLmeanS,biasLstdUS,biasLstdDS]=dB_meanStd(slowTable.BiasLiebe);
    
    [biasITUmeanF,biasITUstdUF,biasITUstdDF]=dB_meanStd(fastTable.BiasITU);
    [biasITUmeanS,biasITUstdUS,biasITUstdDS]=dB_meanStd(slowTable.BiasITU);
    
    [biasLmeanAll,biasLstdUAll,biasLstdDAll]=dB_meanStd(cat(1,fastTable.BiasLiebe,slowTable.BiasLiebe));
    [biasITUmeanAll,biasITUstdUAll,biasITUstdDAll]=dB_meanStd(cat(1,fastTable.BiasITU,slowTable.BiasITU));
    
    disp(['Fast Liebe MM: ',num2str(biasLmeanF),' + ',num2str(biasLstdUF),' - ',num2str(biasLstdDF)]);
    disp(['Fast ITU MM: ',num2str(biasITUmeanF),' + ',num2str(biasITUstdUF),' - ',num2str(biasITUstdDF)]);
    
    disp(['Slow Liebe MM: ',num2str(biasLmeanS),' + ',num2str(biasLstdUS),' - ',num2str(biasLstdDS)]);
    disp(['Slow ITU MM: ',num2str(biasITUmeanS),' + ',num2str(biasITUstdUS),' - ',num2str(biasITUstdDS)]);
    
    disp(['All Liebe MM: ',num2str(biasLmeanAll),' + ',num2str(biasLstdUAll),' - ',num2str(biasLstdDAll)]);
    disp(['All ITU MM: ',num2str(biasITUmeanAll),' + ',num2str(biasITUstdUAll),' - ',num2str(biasITUstdDAll)]);
    
     %% Make plot with overall mean and std
    
    totElevF=[];
    totSig0noBiasF=[];
    
    for ii=1:size(PLTfast,1)
        PLT=PLTfast{ii};
        if isnan(PLT.sig0)
            PLT.sig0noBias=nan(size(PLT.elev));
        end
        totElevF=cat(1,totElevF,PLT.elev);
        totSig0noBiasF=cat(1,totSig0noBiasF,PLT.sig0noBias);
    end
    
    for ii=1:size(PLTslow,1)
        PLT=PLTslow{ii};
        if isnan(PLT.sig0)
            PLT.sig0noBias=nan(size(PLT.elev));
        end
        totElevF=cat(1,totElevF,PLT.elev);
        totSig0noBiasF=cat(1,totSig0noBiasF,PLT.sig0noBias);
    end
    
        elSig=cat(2,totElevF,totSig0noBiasF);
        binInds=discretize(totElevF,edgesFast);
        
        meanSig0noBias=nan(length(edgesFast-1),2);
        stdUSig0noBias=nan(length(edgesFast-1),1);
        stdDSig0noBias=nan(length(edgesFast-1),1);
        
        for jj=1:length(edgesFast-1)
            rowInds=find(binInds==jj);
            data=elSig(rowInds,:);
            meanSig0noBias(jj,1)=nanmean(data(:,1));
            [meanSig0noBias(jj,2),stdUSig0noBias(jj),stdDSig0noBias(jj)]=dB_meanStd(data(:,2));
        end
        
    close all
    
    %fast scans
    
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
     
    errorbar(meanSig0noBias(:,1),meanSig0noBias(:,2),stdDSig0noBias,stdUSig0noBias,'Color','k', 'LineStyle', 'none');
    plot(meanSig0noBias(:,1),meanSig0noBias(:,2),'Color','r','linewidth',2);
    grid on
    
    set(gca,'XLim',[0 25]);
    set(gca,'YLim',[-10,15]);
    xlabel('Incidence Angle, degrees');
    ylabel('Mean raw normalized Radar Cross Section, dB');
    title(['Mean non-corrected Sigma0 vs Incidence Angle, smoothed']);
    set(gcf,'PaperPositionMode','auto')
    print([figdir 'SOCRATES_sigma0noBias_totMean'],'-dpng','-r0')   