% Compare data from before and after velocity correction

clear all
close all

savefig=1;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

plotLine=1;
plotSurf=1;

testCases=readtable(['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/velCorr/inFiles/testCases_cset.dat']);

figdir=['/h/eol/romatsch/hcrCalib/velCorr/',project,'/velFigs_final/'];
formatOut = 'yyyymmdd_HHMM';

indir=HCRdir(project,quality,freqData);

polyTimePeriod=15; %Time period for poly fit in seconds
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
        data.TOPO=[];
%        data.northward_velocity=[];
%        data.eastward_velocity=[];
%        data.azimuth=[];
        
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
            
        %% Fill in extinct echo
        
        % Remove up pointing
        dbzDown=data.DBZ;
        dbzDown(:,data.elevation>0)=nan;
        
        velDown=data.VEL;
        velDown(:,data.elevation>0)=nan;
        
        % Get surface indices
        [linInd rowInd rangeToSurf] = hcrSurfInds(data);
        
        surfDBZ=dbzDown(linInd);
        surfDBZlin=10.^(surfDBZ./10);
        
        % SurfVel is the surface velocity that is used for the polinomial
        % fit
        surfVelOrig=velDown(linInd)';
        surfVel=velDown(linInd);
        
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
        
        aboveGround=data.altitude-data.TOPO;
        
        surfNan=find(isnan(surfVel) & data.elevation<=0 & aboveGround>10);
        
        for ll=1:length(surfNan)
            if ~isnan(surfMean(surfNan(ll)-1)) % At the beginning of the gap we have moving average data
                surfVel(surfNan(ll))=surfMean(surfNan(ll)-1);
            else % Once the moving average turns nan we just keep the previous one going
                surfVel(surfNan(ll))=surfVel(surfNan(ll)-1);
            end
        end
        %% Make poly fit
        velSmoothPorig=vel2vel_corr_testPoly(surfVel,data.time,polyTimePeriod,polyOrder);
        
        % Fill in data where we don't have polyfit data
        velSmoothP=movmedian(surfVelOrig,100,'omitnan');
        velSmoothP(~isnan(velSmoothPorig))=velSmoothPorig(~isnan(velSmoothPorig));
        
        % Vel corr
        velCorrP=data.VEL-velSmoothP';
        velCorrSurfP=velCorrP(linInd);
        
        %% Make output
        velCorrOut=data.VEL;
        velCorrOut(:,data.elevation<=0)=velCorrP(:,data.elevation<=0);
        
        %% Line plot
        
        close all
        
        if plotLine
            f1=figure('DefaultAxesFontSize',12);
            set(f1,'Position',[200 500 1500 1000]);
            set(f1,'renderer','painters');
            
            subplot(2,1,1)
            
            colorIn=lines(length(polyTimePeriod));
            legIn={'VelSurfAngCorr','VelSurfAngCorrFilled','VelSurfPolyFit'};
            
            hold on
            plot(data.time,data.VEL(linInd),'color',[0.5 0.5 0.5],'linewidth',1.5);
            plot(data.time,surfVel,'-c','linewidth',1.5);
            plot(data.time,velSmoothP,'-b','linewidth',1.5);
            ylim([nanmedian(data.VEL(linInd))-1 nanmedian(data.VEL(linInd))+1]);
            ylabel('Velocity [m/s]');
            legend(legIn,'Orientation','horizontal');
            xlim([startTime,endTime]);
            ax=gca;
            ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
            grid on
            grid minor
            
            subplot(2,1,2)
            
            legIn2=[];
            
            hold on
            plot(data.time,velCorrSurfP,'g','linewidth',1);
            
            legend('Poly Corr','Orientation','horizontal');
            xlim([startTime,endTime]);
            ax=gca;
            ax.XAxis.MinorTickValues = [startTime:seconds(10):endTime];
            grid on
            grid minor
            ylabel('Velocity [m/s]');
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f1, [figdir,'lines_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
            end
        end
        %% Plot vel field
        
        if plotSurf
            
            f2=figure('DefaultAxesFontSize',12);
            set(f2,'Position',[200 500 1500 1200]);
            
            % Vel
            subplot(2,1,1)
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
            
            title('Vel');
            
            % VelCorr
            
            subplot(2,1,2)
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
            
            title('VelCorr');
            
            if savefig
                set(gcf,'PaperPositionMode','auto')
                print(f2, [figdir,'velAll_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
            end
        end
    end
end