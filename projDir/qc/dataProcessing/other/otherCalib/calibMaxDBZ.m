% Calibrate radar using bright band method

clear all;
close all;

%caseIn=11;

addpath('/h/eol/romatsch/git/private/utils/')

filedir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/';

figdir='/h/eol/romatsch/hcrCalib/otherCalib/figs/';

plotSingleFigs=1;

% The input file has columns
%startTime (yyyy mm dd HH MM SS)
%endTime
%minimum altitude
%maximum altitu

types='Brightband';
%types='Liquid';
%types='IcephaseUp';

inlist=readtable(['/h/eol/romatsch/hcrCalib/otherCalib/inFiles/caseFile',types,'.dat']);

maxReflAll=[];

for caseIn=1:size(inlist,1)
    
    disp(['Case ',num2str(caseIn),' from ',num2str(size(inlist,1))]);
    
    startTime=datetime(inlist{caseIn,1:6});
    endTime=datetime(inlist{caseIn,7:12});
    
    altMin=inlist{caseIn,13};
    altMax=inlist{caseIn,14};
        
    fileList=makeFileList(filedir,startTime,endTime,1);
    
    outstring1=datestr(startTime,'yyyymmdd_HHMMSS');
    outstring2=datestr(endTime,'yyyymmdd_HHMMSS');
    
    RADIUS  = 6378.140; % radius of earth, in km
    RC = 4*RADIUS; % radius of curvature for "4/3 earth" model
    
    PLT = struct();
    PLT.time = [];
    PLT.refl = [];  % max reflectivity
    PLT.vel  = [];
    PLT.lon = [];
    PLT.lat = [];
    PLT.alt = [];
    
    for ii=1:size(fileList,2)
        indata=fileList{ii};
        
        startTimeIn=ncread(indata,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(indata,'time')';
        
        timeIn=startTimeFile+seconds(timeRead);
        
        %calculate altitude data
        rangeIn=ncread(indata,'range');
        azimuthIn=ncread(indata,'azimuth');
        elevIn=ncread(indata,'elevation');
        latIn=ncread(indata,'latitude');
        lonIn=ncread(indata,'longitude');
        altIn=ncread(indata,'altitude');
        
        if length(elevIn)==1 && length(azimuthIn)>1
            elevIn = repmat(elevIn,size(azimuthIn));
        elseif length(azimuthIn)==1 && length(elevIn)>1
            azimuthIn = repmat(azimuthIn,size(elevIn));
        elseif length(azimuthIn)~=length(elevIn)
            error('length of theta should equal length of elev')
        end
        [thetas,rs] = ndgrid(azimuthIn,rangeIn);
        [elevs,rs] = ndgrid(elevIn,rangeIn);
        % needed to put 0 for OALT since it would be counted twice
        
        elev = mod(elevs,360); % turn elevs to be between 0 and 359.999
        flip_elev_inds = abs(elev-180)<90;
        elev(flip_elev_inds) = 180-elev(flip_elev_inds);
        
        % convert az, elev into radians
        az = thetas * pi/180;
        elev = elev * pi/180;
        
        % compute new coordinates
        xydist = RC*(sin(rs/RC-elev) + sin(elev));
        x = sin(az).*xydist;
        y = cos(az).*xydist;
        z = RC*(cos(rs/RC-elev) - cos(elev));
        
        %  flip_elev_inds = resize(flip_elev_inds,size(x));
        x(flip_elev_inds) = -1*x(flip_elev_inds);
        y(flip_elev_inds) = -1*y(flip_elev_inds);
        
        olat=latIn;
        olon=lonIn;
        oalt=altIn;
        
        % convert olat, olon into radians
        olat = olat * pi/180;
        olon = olon * pi/180;
        
        % compute new coordinates
        lats = asin( (y.*cos(olat) + (RADIUS+oalt+z).*sin(olat)) ./ ...
            sqrt(x.^2+y.^2+(RADIUS+oalt+z).^2)) * 180/pi;
        lons = atan2( (x.*cos(olon) + sin(olon).*((RADIUS+oalt+z).*cos(olat) - y.*sin(olat))), ...
            (-x.*sin(olon) + cos(olon).*((RADIUS+oalt+z).*cos(olat) - y.*sin(olat))) ) * 180/pi;
        alts = sqrt(x.^2 + y.^2 + (RADIUS+oalt+z).^2) - RADIUS;
        
        PLT.vel=cat(2, PLT.vel, ncread(indata,'VEL'));
        PLT.refl=cat(2, PLT.refl, ncread(indata,'DBZ'));
        PLT.alt   = cat(2, PLT.alt,  alts');
        PLT.lat   = cat(2,PLT.lat, lats');
        PLT.lon   = cat(2,PLT.lon,  lons');
        PLT.time = [ PLT.time;  timeIn' ];
    end
    
    
    %% Bright band or ice phase plot
    close all
    
    if ~isnan(altMin)
        
        % Get data within box
        timeInds=find(PLT.time>=startTime & PLT.time<=endTime);
        
        timeBox=PLT.time(timeInds);
        reflTime=PLT.refl(:,timeInds);
        altTime=PLT.alt(:,timeInds);
        
        altInds=find(altTime<altMin | altTime>altMax);
        reflBox=reflTime;
        reflBox(altInds)=nan;
        
        [maxRefl maxGate]=nanmax(reflBox,[],1);
        
        maxReflAll=cat(2,maxReflAll,maxRefl);
        
        if plotSingleFigs
        
        figure;
        
        subplot(3,1,1)
        fig1=surf(PLT.time,PLT.alt,PLT.refl);
        fig1.EdgeColor='none';
        view(2);
        set(gcf,'Position',[200 400 1400 1200]);
        
        fld=fig1.CData;
        
        col_def1 = nan(size(fld));
        col_def2 = nan(size(fld));
        col_def3 = nan(size(fld));
        
        limits=[-inf (-43:3:23) inf];
        color_map=colormap(dbz_default);
        
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
        fig1.CData=col_def;
        
        colormap(gca,color_map);
        caxis([0 size(color_map,1)]);
        caxis_yticks=(1:1:size(color_map,1)-1);
        caxis_ytick_labels=num2str(limits(2:end-1)');
        while length(caxis_yticks)>16
            caxis_yticks=caxis_yticks(1:2:end);
            caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
        end
        
        ylim([-1000 8000]);
        ylabel('Altitude (m)');
        xlim([startTime endTime]);
        
        % Find bright band
        hold on;
        
        zplot = get(fig1,'ZData');
        set(fig1,'ZData',zplot-10)
        % z_max = max(max(get(fig1,'Zdata')))
        % line(1:50,1:50,z_max*ones(1,50))
        
        line([PLT.time(timeInds(1)),PLT.time(timeInds(end))],[altMin,altMin],[1000 1000],'color','c','linewidth',2);
        line([PLT.time(timeInds(1)),PLT.time(timeInds(end))],[altMax,altMax],[1000 1000],'color','c','linewidth',2);
        line([PLT.time(timeInds(1)),PLT.time(timeInds(1))],[altMin,altMax],[1000 1000],'color','c','linewidth',2);
        line([PLT.time(timeInds(end)),PLT.time(timeInds(end))],[altMin,altMax],[1000 1000],'color','c','linewidth',2);
        
        subplot(3,1,3)
        plot(timeBox,maxRefl);
        grid on;
        ylabel('Max. refl. (dBZ)');
        ylim([0 20]);
        xlim([startTime endTime]);
        
        subplot(3,1,2)
        fig2=surf(PLT.time,PLT.alt,PLT.refl);
        fig2.EdgeColor='none';
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
        set(get(hcb,'Title'),'String','dBZ');
        set(hcb,'location','northoutside');
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
        
        ylim([altMin altMax]);
        ylabel('Altitude (m)');
        xlim([startTime endTime]);
        
        set(gcf,'PaperPositionMode','auto')
        print([figdir types '_' outstring1 '_to_' outstring2],'-dpng','-r0')
    end
    end
    
end

%% Histogram

above15=maxReflAll;
above15(find(above15<15))=[];

f1=figure;
set(f1,'Position',[200 100 800 800],'DefaultAxesFontSize',14);

edges=(15:0.5:24);
histogram(above15,edges);

xlim([15 23.5]);
ylim([0 1500]);

xlabel('Max. refl. (dBZ)');
ylabel('Number of samples');
title([types ' histogram of reflectivities > 15 dBZ'],'interpreter','none');

set(gcf,'PaperPositionMode','auto')
print([figdir types '_histogram'],'-dpng','-r0')