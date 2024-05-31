function data=interp_HCR2grid(dataRadial,fig);
% Interpolate HCR data to rectangular grid
data=[];

dataVars=fieldnames(dataRadial);

% Set up grid
maxAsl=max(max(dataRadial.asl));
newAsl=0:dataRadial.range(2)-dataRadial.range(1):maxAsl;
newAsl=double(newAsl);

aslMat=repmat(double(newAsl'),1,length(dataRadial.time));
timeMat=repmat(datenum(dataRadial.time),size(aslMat,1),1);

data.asl=aslMat;
data.range=dataRadial.range(:,1);

for ii=1:length(dataVars);
    dataIn=dataRadial.(dataVars{ii});
    
    % 2D field
    if min(size(dataIn))>1 & ~strcmp(dataVars{ii},'range') & ~strcmp(dataVars{ii},'asl')
        dataInCensored=dataIn;
        dataInCensored(isnan(dataIn))=-99999;
        vq=[];
        for mm=1:size(timeMat,2)
            x=dataRadial.asl(:,mm);
            y=dataInCensored(:,mm);
            xy=cat(2,x,y);
            nanXY=find(any(isnan(xy),2));
            xy(nanXY,:)=[];
            if ~isempty(xy)
                vq1 = interp1(xy(:,1),xy(:,2),newAsl,'nearest');
                vq=cat(2,vq,vq1');
            else
                vq=cat(2,vq,nan(length(newAsl),1));
            end
        end
        vq(vq==-99999)=nan;
        data.(dataVars{ii})=vq;
        
        if fig
            close all
            fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
            colormap jet;
            
            ax1=subplot(2,1,1);
            hold on;
            
            sub1=surf(dataRadial.time,dataRadial.asl./1000,dataIn,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            xlim([dataRadial.time(1),dataRadial.time(end)]);
            ylim([-0.2 15]);
            title([dataVars{ii},' radial']);
            grid on
            
            ax2=subplot(2,1,2);
            
            hold on;
            sub2=surf(datetime(timeMat,'ConvertFrom','datenum'),aslMat./1000,vq,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            xlim([dataRadial.time(1),dataRadial.time(end)]);
            ylim([-0.2 15]);
            title([dataVars{ii},' gridded']);
            grid on
        end
    elseif ~strcmp(dataVars{ii},'range') & ~strcmp(dataVars{ii},'asl') % 1D field -> just copy data
        data.(dataVars{ii})=dataRadial.(dataVars{ii});
    end
end
end

