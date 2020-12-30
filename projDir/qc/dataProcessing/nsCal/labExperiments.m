% Perform noise cal

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/labCal/'];

indir='~/git/HCR_configuration/projDir/calibration/data/lna_temp/';

inlist={'TsPrint.lna_temp.20201105.txt';
    'TsPrint.lna_temp.20201105.txt';
    'TsPrint.lna_temp.20201119.txt';
    'TsPrint.lna_temp.20201216.txt';
    %'TsPrint.lna_temp.20201216.txt';
    'TsPrint.lna_temp.20201217.txt';
    'TsPrint.lna_temp.20201218.txt'};

fileTimes=[datetime(2020,11,5,20,25,0),datetime(2020,11,5,20,31,0);
    datetime(2020,11,5,20,49,0),datetime(2020,11,5,20,55,0);
    datetime(2020,11,19,21,3,0),datetime(2020,11,19,21,8,0);
    datetime(2020,12,16,22,26,0),datetime(2020,12,16,22,30,0);
    %datetime(2020,12,16,22,40,0),datetime(2020,12,16,22,45,0);
    datetime(2020,12,17,16,17,0),datetime(2020,12,17,16,21,0);
    datetime(2020,12,18,22,12,0),datetime(2020,12,18,22,30,0)];

timeLag=[7,9;
    7,9;
    5,9;
    5,9;
    %9,10;
    9,9;
    9,8];

%% Loop through cases

for ii=1:size(inlist,1)
    
    disp(datestr(fileTimes(ii,1),'yyyy-mm-dd HH:MM:SS'));
    
    %% Read data
    fid=fopen([indir,inlist{ii}],'r');
    slurp=fscanf(fid,'%c');
    fclose(fid);
    
    M=strread(slurp,'%s','delimiter','\n');
    
    commentInds = find(contains(M,'#'));
    M(commentInds)=[];
    
    time=NaT(length(M),1);
    hc=nan(length(M),1);
    vc=nan(length(M),1);
    ht=nan(length(M),1);
    vt=nan(length(M),1);
    
    for jj=1:length(M)
        lineIn=M{jj};
        lsplit=split(lineIn,' ');
        lsplit=lsplit(~cellfun('isempty',lsplit));
        
        timeCell=lsplit{1};
        time(jj)=datetime(str2num(timeCell(1:4)),str2num(timeCell(6:7)),str2num(timeCell(9:10)),...
            str2num(timeCell(12:13)),str2num(timeCell(15:16)),str2num(timeCell(18:end)));
        
        hc(jj)=str2num(lsplit{5});
        vc(jj)=str2num(lsplit{8});
        ht(jj)=str2num(lsplit{13});
        vt(jj)=str2num(lsplit{14});
    end
    
    %% Clean up data
    % Only use 25-35 C temperatures       
    tinds=find(time>=fileTimes(ii,1) & time<=fileTimes(ii,2));
    ht=ht(tinds);
    hc=hc(tinds);
    htime=time(tinds);
    
    vt=vt(tinds);
    vc=vc(tinds);
    vtime=time(tinds);
    
    %% Adjust times
    TTHC=timetable(htime,hc);
    TTHT=timetable(htime,ht);
    TTHCadj=timetable(htime+seconds(timeLag(ii,1)),hc);
        
    TTVC=timetable(vtime,vc);
    TTVT=timetable(vtime,vt);
    TTVCadj=timetable(vtime+seconds(timeLag(ii,2)),vc);
        
    TTHall=synchronize(TTHT,TTHCadj,'first','nearest');    
    TTHall(find(TTHall.htime<TTHCadj.Time(1)),:)=[];
    
    TTVall=synchronize(TTVT,TTVCadj,'first','nearest');    
    TTVall(find(TTVall.vtime<TTVCadj.Time(1)),:)=[];
    
    % Plot
    close all
   
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,820]);
    
    subplot(2,1,1)
    hold on
    plot(htime,hc,'-m','linewidth',1);
    plot(TTHall.htime,TTHall.hc,'-r','linewidth',1.5);
    ylabel('Power (dBm)');
    
    yyaxis right
    hold on
    plot(TTHall.htime,TTHall.ht,'-c','linewidth',1);
    plot(TTHall.htime,movmean(TTHall.ht,20),'-b','linewidth',1.5);
    ylabel('Temperature(C)');
    set(gca,'YColor','k');
    
    title(['H power vs temperature ',datestr(TTHall.htime(1),'yyyy-mm-dd HH:MM:SS')]);
    grid on
    
    legend({'H power','H power adj','H temp raw','H temp smooth'},'location','southeast');
    xlim([TTHall.htime(1),TTHall.htime(end)]);
    
    subplot(2,1,2)
    hold on
    plot(vtime,vc,'-m','linewidth',1);
    plot(TTVall.vtime,TTVall.vc,'-r','linewidth',1.5);
    ylabel('Power (dBm)');
    
    yyaxis right
    hold on
    plot(TTVall.vtime,TTVall.vt,'-c','linewidth',1);
    plot(TTVall.vtime,movmean(TTVall.vt,20),'-b','linewidth',1.5);
    ylabel('Temperature(C)');
    set(gca,'YColor','k');
    
    title(['V power vs temperature ',datestr(TTVall.vtime(1),'yyyy-mm-dd HH:MM:SS')]);
    grid on
    
    legend({'V power','V power adj','V temp raw','V temp smooth'},'location','southeast');
    xlim([TTVall.vtime(1),TTVall.vtime(end)]);
    
    print([figdir,'powerVStemp_',datestr(TTVall.vtime(1),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    %% Fit
   
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,600,920]);
    
    subplot(2,1,1)
    hold on
    scatter(TTHall.ht,TTHall.hc,'filled');
    xlimits=xlim;
    ylimits=ylim;
    
    fitOrth=gmregress(TTHall.ht,TTHall.hc,1);
    fitAll=[fitOrth(2) fitOrth(1)];
    xFit = xlimits(1):0.1:xlimits(2);
    yFit = polyval(fitAll, xFit);
    
    plot(xFit, yFit,'-r','linewidth',2);
    
    xlim(xlimits);
    ylim(ylimits);
    grid on
    
    xlabel('Temperature [C]');
    ylabel('H power [dBm]');
        
    text(xlimits(1)+0.5,ylimits(2)-(ylimits(2)-ylimits(1))/10,...
        ['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
    
    title(['H power vs temperature ',datestr(TTHall.htime(1),'yyyy-mm-dd HH:MM:SS')]);
    
    subplot(2,1,2)
    hold on
    scatter(TTVall.vt,TTVall.vc,'filled');
    xlimits=xlim;
    ylimits=ylim;
    
    fitOrth=gmregress(TTVall.vt,TTVall.vc,1);
    fitAll=[fitOrth(2) fitOrth(1)];
    xFit = xlimits(1):0.1:xlimits(2);
    yFit = polyval(fitAll, xFit);
    
    plot(xFit, yFit,'-r','linewidth',2);
    
    xlim(xlimits);
    ylim(ylimits);
    grid on
    
    xlabel('Temperature [C]');
    ylabel('V power [dBm]');
        
    text(xlimits(1)+0.5,ylimits(2)-(ylimits(2)-ylimits(1))/10,...
        ['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
    
    title(['V power vs temperature ',datestr(TTVall.vtime(1),'yyyy-mm-dd HH:MM:SS')]);
    
    print([figdir 'scatter_',datestr(TTVall.vtime(1),'yyyymmdd_HHMMSS')],'-dpng','-r0')
end
