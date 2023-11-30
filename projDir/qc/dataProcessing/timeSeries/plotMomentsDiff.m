function plotMomentsDiff(momentsTime,momentsSpec,figdir,project,ylimUpper,showPlot)
momentsTime=rmfield(momentsTime,'azimuth_vc');
momentsTime=rmfield(momentsTime,'eastward_velocity');
momentsTime=rmfield(momentsTime,'elevation');
momentsTime=rmfield(momentsTime,'northward_velocity');
momentsTime=rmfield(momentsTime,'range');
momentsTime=rmfield(momentsTime,'time');
momentsTime=rmfield(momentsTime,'vertical_velocity');

tFields=fields(momentsTime);
sFields=fields(momentsSpec);

commFields=intersect(tFields,sFields);

for ii=1:length(commFields)

    if ~isempty(find(~isnan(momentsSpec.(commFields{ii}))))

        f1 = figure('Position',[200 500 900 900],'DefaultAxesFontSize',12,'visible',showPlot);

        colormap jet
        s1=subplot(3,1,1);

        hold on
        surf(momentsSpec.time,momentsTime.asl./1000,momentsTime.(commFields{ii}),'edgecolor','none');
        view(2);
        ylabel('Range (km)');
        ylim([0 ylimUpper]);
        xlim([momentsSpec.time(1),momentsSpec.time(end)]);
        colorbar
        grid on
        box on
        title([commFields{ii},' time domain']);

        s2=subplot(3,1,2);

        hold on
        surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.(commFields{ii}),'edgecolor','none');
        view(2);
        ylabel('Range (km)');
        ylim([0 ylimUpper]);
        xlim([momentsSpec.time(1),momentsSpec.time(end)]);
        colorbar
        grid on
        box on
        title([commFields{ii},' spectral domain']);

        s3=subplot(3,1,3);

        hold on
        surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.(commFields{ii})-momentsTime.(commFields{ii}),'edgecolor','none');
        view(2);
        ylabel('Range (km)');
        clim([-0.2 0.2]);
        ylim([0 ylimUpper]);
        xlim([momentsSpec.time(1),momentsSpec.time(end)]);
        colorbar
        grid on
        box on
        title([commFields{ii},' spectral minus time domain']);

        if strcmp(commFields{ii},'dbz')
            s1.CLim=[-60,20];
            s2.CLim=[-60,20];
        elseif strcmp(commFields{ii},'velRaw') | strcmp(commFields{ii},'vel')
            s1.CLim=[-8,8];
            s2.CLim=[-8,8];
        elseif strcmp(commFields{ii},'powerV')
            s1.CLim=[-110,-40];
            s2.CLim=[-110,-40];
        elseif strcmp(commFields{ii},'powerH')
            s1.CLim=[-110,-40];
            s2.CLim=[-110,-40];
        elseif strcmp(commFields{ii},'width')
            s1.CLim=[0,4];
            s2.CLim=[0,4];
        elseif strcmp(commFields{ii},'snr')
            s1.CLim=[-20,70];
            s2.CLim=[-20,70];
        elseif strcmp(commFields{ii},'ldr')
            s1.CLim=[-40,10];
            s2.CLim=[-40,10];
        elseif strcmp(commFields{ii},'skew')
            s1.CLim=[-1,1];
            s2.CLim=[-1,1];
        elseif strcmp(commFields{ii},'kurt')
            s1.CLim=[0,10];
            s2.CLim=[0,10];
        end

        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_diff_',(commFields{ii}),'_',datestr(momentsSpec.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momentsSpec.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    end
end
end