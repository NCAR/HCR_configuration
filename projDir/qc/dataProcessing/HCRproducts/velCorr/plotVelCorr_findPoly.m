% Compare data from before and after velocity correction

clear all
close all

savefig=1;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

testCases=readtable(['~/git/HCR_configuration/projDir/qc/dataProcessing/velCorr/inFiles/testCases_',project,'.dat']);

figdir=['/h/eol/romatsch/hcrCalib/velCorr/',project,'/velFigs/'];
formatOut = 'yyyymmdd_HHMM';

indir=HCRdir(project,quality,freqData);

polyTimePeriod=[10 15 20]; %Vector of time periods for poly fit in seconds
polyUsed=2; % Index used for correction from polyTimePeriod
polyOrder=3; % Order of polynomial fit

maxEdge=14.6; % Upper edge for plotting
color_map=colormap(vel_default(29));

limits=-4.05:0.3:4.05;
limits=[-inf limits inf];

for kk=1:size(testCases,1)
    startTime=datetime(testCases{kk,1},testCases{kk,2},testCases{kk,3}, ...
        testCases{kk,4},testCases{kk,5},0);
    endTime=datetime(testCases{kk,6},testCases{kk,7},testCases{kk,8}, ...
        testCases{kk,9},testCases{kk,10},0);
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        data=[];
        
        data.vertical_velocity=[];
        data.nyquist_velocity=[];
        data.DBZ=[];
        data.VEL=[];
        data.VEL_RAW=[];
        data.TOPO=[];
        data.northward_velocity=[];
        data.eastward_velocity=[];
        data.azimuth=[];
        
        dataVars=fieldnames(data);
        
        if length(fileList)==0
            disp('No data files found.');
            return
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
        
        dataVars=dataVars(~cellfun('isempty',dataVars));
        
        %% Correct angles
        
        % Get surface indices
        [linInd rowInd rangeToSurf] = hcrSurfInds(data);
        
%         % Use these equations when folding is a problem (e.g. in ocean scans)
%         % For VEL
%         posInds=find(data.VEL(linInd)>0);
%         data.VEL(:,posInds)=data.VEL(:,posInds)-2.*data.nyquist_velocity(posInds);
% %         % For VEL_RAW
% %         posInds2=find(data.VEL_RAW(linInd)>0);
% %         data.VEL_RAW(:,posInds2)=data.VEL_RAW(:,posInds2)-2.*data.nyquist_velocity(posInds2);
        
%         xxCorr=sind(data.azimuth).*cosd(data.elevation).*data.eastward_velocity;
%         yyCorr=cosd(data.azimuth).*cosd(data.elevation).*data.northward_velocity;
%         
%         % Use this equations when starting from VEL
%         velAngCorr=data.VEL+xxCorr+yyCorr;
%         
%         % Use these equations when starting from VEL_RAW
%         zzCorr=sind(data.elevation).*data.vertical_velocity;
%         velAngCorr2=data.VEL_RAW+xxCorr+yyCorr+zzCorr;
            
        velAngCorr=data.VEL;
        %% Fill in extinct echo
        surfDBZ=data.DBZ(linInd);
        surfDBZlin=10.^(surfDBZ./10);
        
        % SurfVel is the surface velocity that is used for the polinomial
        % fit
        surfVel=velAngCorr(linInd);
        
        % We remove all the data that we don't want to include in the fit
        
        % Remove all data where we don't see the surface
        surfVel(isnan(surfDBZlin))=nan;
        
        % Remove data that is out of range
        surfVel(isnan(rowInd))=nan;
        
        % Calculate standard deviation
        surfStd=movstd(surfVel,100,'includenan'); % Sets to nan when there is at least one nan
        
        % Enlarge missing data areas
        surfTemp=movmean(surfVel,100,'includenan'); % Sets to nan when there is at least one nan
        surfVel(isnan(surfTemp))=nan;
        
        % Remove data points where the surface echo is weak and where the
        % standard deviation is too high
        surfVel(surfDBZlin<10000 & surfStd>0.5)=nan;
                      
        % Fill in the missing data with moving average from before the gap
        surfMean=movmedian(surfVel,100,'omitnan'); % Ignores nans and uses the rest of the data
        surfMean(isnan(surfVel))=nan;
        
        surfNan=find(isnan(surfVel));
        
        for ll=1:length(surfNan)
            if ll==1
                surfVel(surfNan(ll))=surfMean(surfNan(ll));
            elseif ~isnan(surfMean(surfNan(ll)-1)) % At the beginning of the gap we have moving average data
                surfVel(surfNan(ll))=surfMean(surfNan(ll)-1);
            else % Once the moving average turns nan we just keep the previous one going
                surfVel(surfNan(ll))=surfVel(surfNan(ll)-1);
            end
        end
        %% Make poly fit
        velSmoothP=nan(length(polyTimePeriod),length(linInd));
        for ii=1:length(polyTimePeriod)
            velSmoothP(ii,:)=vel2vel_corr_testPoly(surfVel,data.time,polyTimePeriod(ii),polyOrder);
        end
               
        %% Vel corr
       
        velCorrP=nan(size(data.VEL,1),size(data.VEL,2),length(polyTimePeriod));
        velCorrSurfP=nan(size(data.VEL,2),length(polyTimePeriod));
        for jj=1:length(polyTimePeriod)
            velCorrT=velAngCorr-velSmoothP(jj,:);
            velCorrP(:,:,jj)=velCorrT;
            velCorrSurfP(:,jj)=velCorrT(linInd);
        end
        
        %% Put angle corrected data back in when pointing up
        velCorrOut=velCorrP(:,:,polyUsed);
        velCorrOut(:,data.elevation>0)=velAngCorr(:,data.elevation>0);
        
        %% Line plot
        
        close all
        f1=figure('DefaultAxesFontSize',12);
        set(f1,'Position',[200 500 1500 1500]);
        set(f1,'renderer','painters');
        
        subplot(3,1,1)
        hold on
        plot(data.time,data.VEL_RAW(linInd),'-g','linewidth',1.5);
        plot(data.time,data.VEL(linInd),'color',[0.5 0.5 0.5],'linewidth',1.5);
        ylim([nanmedian(data.VEL(linInd))-1.5 nanmedian(data.VEL(linInd))+1.5]);
        ylabel('Velocity [m/s]');
        legend('VelRawSurf','VelSurf','Orientation','horizontal');
        title([datestr(startTime,formatOut) ' to ' datestr(endTime,formatOut)],'interpreter','none');
        xlim([startTime,endTime]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        
        subplot(3,1,2)
        
        colorIn=lines(length(polyTimePeriod));
        legIn={'VelSurf','VelSurfAngCorrFilled'};
        
        hold on
        plot(data.time,velAngCorr(linInd),'color',[0.5 0.5 0.5],'linewidth',1.5);
        plot(data.time,surfVel,'-c','linewidth',1.5);
        for ii=1:length(polyTimePeriod)
            plot(data.time,velSmoothP(ii,:),'color',colorIn(ii,:),'linewidth',1.5);
            legIn=[legIn,num2str(polyTimePeriod(ii))];
        end
        ylim([nanmedian(data.VEL(linInd))-1 nanmedian(data.VEL(linInd))+1]);
        ylabel('Velocity [m/s]');
        legend(legIn,'Orientation','horizontal');
        xlim([startTime,endTime]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        
        subplot(3,1,3)
        
        legIn2=[];
        
        hold on
        for ii=1:length(polyTimePeriod)
            plot(data.time,velCorrSurfP(:,ii),'color',colorIn(ii,:),'linewidth',1);
            legIn2{end+1}=num2str(polyTimePeriod(ii));
        end
        legend(legIn2,'Orientation','horizontal');
        xlim([startTime,endTime]);
        ax=gca;
        ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
        grid on
        grid minor
        ylabel('Velocity [m/s]');
        title('Vel Surf Poly Corrected');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f1, [figdir,'lines_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
        end
        
        %% Plot vel field
               
        f2=figure('DefaultAxesFontSize',12);
        set(f2,'Position',[200 500 1500 1500]);
        
        % Vel
        subplot(3,1,1)
        hold on
        fig2=surf(data.time,data.asl./1000,data.VEL_RAW);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
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
        set(get(hcb,'Title'),'String','m/s');
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
        
        title('VEL_RAW','interpreter','none');
        
        % Vel angle corrected
        
        subplot(3,1,2)
        hold on
        fig2=surf(data.time,data.asl./1000,data.VEL);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
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
        set(get(hcb,'Title'),'String','m/s');
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
        
        title('VEL');
        
        % VelCorr
        
        subplot(3,1,3)
        hold on
        fig2=surf(data.time,data.asl./1000,velCorrOut);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
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
        set(get(hcb,'Title'),'String','m/s');
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
  
        title(['VEL_CORR ',num2str(polyTimePeriod(polyUsed)),' s'],'interpreter','none');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f2, [figdir,'velAll_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
        end
        
        %% Make poly fit plot
               
        f3=figure('DefaultAxesFontSize',12);
        set(f3,'Position',[1700 500 1500 1500]);
        
        for jj=1:length(polyTimePeriod)
            
            subplot(length(polyTimePeriod),1,jj)
            hold on
            fig4=surf(data.time,data.asl./1000,velCorrP(:,:,jj));
            fig4.EdgeColor='none';
            ylim([-0.2 maxEdge]);
            xlim([startTime,endTime]);
            view(2);
            
            fld=fig4.CData;
            
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
            fig4.CData=col_def;
            
            hcb=colorbar;
            set(get(hcb,'Title'),'String','m/s');
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
            
            title(['VelCorr poly fit 3, ',num2str(polyTimePeriod(jj)),' s']);
        end
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f3, [figdir,'velPoly_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
        end
    end
end