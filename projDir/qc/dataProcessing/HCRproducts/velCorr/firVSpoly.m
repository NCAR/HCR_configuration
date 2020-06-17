% Compare data from before and after velocity correction

clear all
close all

savefig=0;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

testCases=readtable(['/h/eol/romatsch/hcrCalib/velCorr/otrec/testCases.dat']);

polyTimePeriod=10;
polyOrder=3;

for kk=1:size(testCases,1)
    startTime=datetime(testCases{kk,1},testCases{kk,2},testCases{kk,3}, ...
        testCases{kk,4},testCases{kk,5},0);
    endTime=datetime(testCases{kk,6},testCases{kk,7},testCases{kk,8}, ...
        testCases{kk,9},testCases{kk,10},0);
         
    indir='/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/';
    figdir=['/h/eol/romatsch/hcrCalib/velCorr/otrec/velFigs/scott/'];
    formatOut = 'yyyymmdd_HHMM';
     
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        data.time=[];
        data.alt=[];
        data.elev=[];
        data.range=[];
        data.dbz=[];
        data.vel=[];
          
        % Get uncorrected data
        for ii=1:size(fileList,2)
            infile=fileList{ii};
            
            startTimeIn=ncread(infile,'time_coverage_start')';
            startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
                str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
            timeRead=ncread(infile,'time')';
            timeIn=startTimeFile+seconds(timeRead);
            data.time=[data.time,timeIn];
            
            rangeIn=ncread(infile,'range');
            rangeMat=repmat(rangeIn,1,length(timeIn));
            data.range=[data.range,rangeMat];
            
            data.elev=[data.elev,ncread(infile,'elevation')'];
            
            data.alt=[data.alt,ncread(infile,'altitude')'];
            
            data.dbz=[data.dbz,ncread(infile,'DBZ')];
            
            velIn=ncread(infile,'VEL');
            data.vel=[data.vel,velIn];
        end
  
        %% Get Scott's data
        [velCorrS,velCorrSurfS,rangeToSurfS,velSmoothS]=vel2vel_corr(data.vel,data.dbz,data.range,data.alt);
        
         %% Poly fit
         
         data.DBZ=data.dbz;
         data.altitude=data.alt;
         data.elevation=data.elev;
         data.VEL=data.vel;
         [linInd rowInd rangeToSurf] = hcrSurfInds(data);
        
        %% Interpolate over extinct echo
        surfDBZ=data.DBZ(linInd);
        surfDBZlin=10.^(surfDBZ./10);
        
        surfVel=data.VEL(linInd);
        surfVel(isnan(surfDBZlin))=nan;
        surfStd=movstd(surfVel,100,'omitnan');
        surfVel(surfDBZlin<10000 & surfStd>0.5)=nan;
        %surfVel(surfDBZlin<10000)=nan;
                
        surfMean=movmedian(surfVel,100,'omitnan');
        surfMean(isnan(surfVel))=nan;
        
        surfNan=find(isnan(surfVel));
        
        for ll=1:length(surfNan)
            if ~isnan(surfMean(surfNan(ll)-1))
                surfVel(surfNan(ll))=surfMean(surfNan(ll)-1);
            else
                surfVel(surfNan(ll))=surfVel(surfNan(ll)-1);
            end
        end
        
        %% Make poly fit
        %velSmoothP=nan(length(linInd),1);
        velSmoothP=vel2vel_corr_testPoly(surfVel,data.time,polyTimePeriod,polyOrder);
        
        %% Vel corr
        
        velCorrP=nan(size(data.VEL,1),size(data.VEL,2));
        velCorrSurfP=nan(size(data.VEL,2),1);
        
        velMat=repmat(velSmoothP,1,size(data.VEL,1))';
        velCorrT=data.VEL-velMat';
        velCorrP=velCorrT;
        velCorrSurfP(:,jj)=velCorrT(linInd);
        
        
        %% Plot
        
        close all
        
        % Plot vel field
        
        %calculate above sea level altitudes
        asl=-1*((data.range.*cosd(abs(data.elev)-90)./1000)-data.alt./1000);
        
        maxAlt=max(max(data.alt))/1000;
        maxEdge=ceil( maxAlt/0.5 ) * 0.5;
        
        f2=figure('DefaultAxesFontSize',14);
        set(f2,'Position',[200 500 1500 1500]);
        
        subplot(3,1,1)
        fig2=surf(data.time,asl,data.vel);
        fig2.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig2.CData;
        
        col_def1 = nan(size(fld));
        col_def2 = nan(size(fld));
        col_def3 = nan(size(fld));
        
        color_map=colormap(vel_default(29));
        %color_map=jet(17);
        
        %limits=[-8 -3.8 -3.3 -2.8 -2.3 -1.8 -1.3 -0.8 -0.3 0.3 0.8 1.3 1.8 2.3 2.8 3.3 3.8 8];
        limits=-4.05:0.3:4.05;
        limits=[-inf limits inf];
        
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
        
        subplot(3,1,2)
        fig2=surf(data.time,asl,data.velCorr);
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
        title('VelCorr Mike');
        
        
        subplot(3,1,3)
        fig3=surf(data.time,asl,velCorrS);
        fig3.EdgeColor='none';
        ylim([-0.2 maxEdge]);
        xlim([startTime,endTime]);
        view(2);
        
        fld=fig3.CData;
        
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
        fig3.CData=col_def;
        
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
        title('VelCorr Scott');
        
        if savefig
            set(gcf,'PaperPositionMode','auto')
            print(f2, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_velCompare'],'-dpng','-r0');
        end
        
    end
end