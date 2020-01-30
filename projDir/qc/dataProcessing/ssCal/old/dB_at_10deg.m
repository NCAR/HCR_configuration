%Find the ocean surface reflectivity at 10 degree beam angle


clear all;
close all;

addpath('/h/eol/romatsch/matlab/radar/functions/');
addpath('/h/eol/romatsch/matlab/radar/subCodes/');

%Radar beam angle (degree) away from 180 deg at which reflectivity is wanted
radarAngle=10;
degPlusMinus=0.3; %Tolerance of radar angle in degree

%Time tolerance in seconds
%timePlusMinus=300;
timePlusMinus=6000;

outFile='/h/eol/romatsch/radar/notes/HCR_bias_elevAng.dat';

radarFileDir='/h/eol/romatsch/matlab/radar/fileLists/radarScans/';

allFileList={'sound_20150722_193012.dat',...
'sound_20150724_165608.dat',...
'sound_20150724_230006.dat',...
'sound_20150801_214509.dat',...
'sound_20150812_201151.dat',...
'sound_20150709_222840.dat',...
'sound_20150709_224819.dat',...
'sound_20150803_170007.dat',...
'sound_20150807_193958.dat',...
'sound_20150807_212000.dat'};

%...read file list
%allFileList=cellstr('sound_20150812_165938.dat');
%allFileList=cellstr('sound_20150724_165608.dat');
%allFileList=cellstr('sound_20150724_230006.dat');
%allFileList=cellstr('sound_20150801_214509.dat');
%allFileList=cellstr('sound_20150812_201151.dat');

%allFileList=cellstr('sound_20150709_222840.dat');
%allFileList=cellstr('sound_20150709_224819.dat');
%allFileList=cellstr('sound_20150803_170007.dat');
%allFileList=cellstr('sound_20150807_193958.dat');
%allFileList=cellstr('sound_20150807_212000.dat');

output=zeros(length(allFileList),8);

for j=1:length(allFileList)
    
    radarFile=allFileList{j};
    
    fileList=importdata(strcat(radarFileDir,radarFile));
    
    %Sounding date
    soundDate=radarFile(7:14);
    
    %Sounding time
    soundHH=radarFile(16:17);
    soundMM=radarFile(18:19);
    soundSS=radarFile(20:21);
    
    direc = '/scr/eldora1/rsfdata/cset/hcr/qc/cfradial/velcorr/10hz/';
    %direc = '/scr/eldora1/rsfdata/cset/hcr/cfradial/moments/10hz/';
    %direc = '/scr/eldora1/rsfdata/cset/hcr/cfradial/moments/100hz/';
    set(0,'DefaultTextInterpreter','none');
    
    PLT = struct();
    PLT.time = [];
    PLT.rota = [];
    PLT.refl = [];  % max reflectivity
    PLT.alt  = [];
    PLT.pitch = [];
    PLT.range = [];
    PLT.elevation = [];
    
    for i=1:length(fileList);
        disp(['File ',num2str(i),' of ',num2str(length(fileList))]);
        f_name = strcat(direc,soundDate,'/',cell2mat(fileList(i)));
        
        [ ATT GATT ] = f_get_hdf5_param_atts( f_name(1,:) );
        %Flds = f_determine_radar_fields( ATT );
        % sFlds = select_params(Flds);
        % tFlds = determine_time_depends( SWP );
        
        SWP = f_get_and_scale_hdf5_data(deblank(f_name), ATT,[fields(ATT)' ]);
        
        ngates = size(SWP.DBZ,1);
        nbeams = size(SWP.DBZ,2);
        
        PLT.rota  = [ PLT.rota;  SWP.rotation ];
        PLT.alt   = [ PLT.alt;   SWP.altitude ];
        PLT.pitch = [ PLT.pitch; SWP.pitch ];
        PLT.refl=cat(2,PLT.refl,SWP.DBZ);
        PLT.range=cat(2,PLT.range,SWP.range);
        PLT.elevation  = [ PLT.elevation;  SWP.elevation ];
        
        clear tmpsec;
        tmpsec = f_get_cfradnc_daysecs( SWP );
        PLT.time = [ PLT.time;  tmpsec ];
        if(size(tmpsec,1) ~= nbeams);
            fprintf('Error on file: nbeams not equal to number of time points\n');
        end;
    end;
    
    f=SWP.frequency/10^9;
    
    if(mean(max(PLT.range,[],2)-min(PLT.range,[],2))==0);
        PLT.range=PLT.range(:,1);
    else
        fprintf('Ranges are not all the same.\n');
    end;
    
    % There are high reflectivities in the first gates that need to be removed.
    badRow=8;
    empty=1;
    
    while (badRow<30 && empty>0)
        empty=length(find(~isnan(PLT.refl(badRow,:))));
        badRow=badRow+1;
    end;
    
    disp(['First row with only nans = ',num2str(badRow)]);
    
    PLT.refl(1:badRow-1,:)=nan;
    
    soundTimeLong=str2num(soundHH)*3600+str2num(soundMM)*60+str2num(soundSS);
    
    diffTime=abs(PLT.time-soundTimeLong);
    timeInd=find(diffTime<=timePlusMinus);
    
    tPLT = struct();
    tPLT.time=PLT.time(timeInd,:);
    tPLT.rota=PLT.rota(timeInd,:);
    tPLT.refl=PLT.refl(:,timeInd);
    tPLT.elevation=PLT.elevation(timeInd,:);
    
    %This part is for radar rotation angle
    %angPlus10=find(tPLT.rota>(180+radarAngle-degPlusMinus) & tPLT.rota<(180+radarAngle+degPlusMinus));
    %angMinus10=find(tPLT.rota>(180-radarAngle-degPlusMinus) & tPLT.rota<(180-radarAngle+degPlusMinus));
    %angInd=cat(1,angPlus10,angMinus10);
    %Only left or right side values
    %angInd=angPlus10;
    %angInd=angMinus10;
    
    %This part is for radar elevation angle
    angInd=find(tPLT.elevation>(-90+radarAngle-degPlusMinus) & tPLT.elevation<(-90+radarAngle+degPlusMinus));
    
    angTime=tPLT.time(angInd,:);
    angRefl=tPLT.refl(:,angInd);
    
    numNonNa=sum(~isnan(angRefl),1);
    maxRefl=max(angRefl);
    
    [reflClean, badIndRefl,outlRefl]=f_deleteoutliers(maxRefl,0.05,1);
    badIndNum=find(abs(numNonNa-mode(numNonNa))>4);
    numClean=numNonNa;
    numClean(badIndNum)=nan;
    %[numClean, badIndNum,outlNum]=f_deleteoutliers(numNonNa,0.05,1);
    
    close all;
    
    figure
    plot(numNonNa,'b');
    hold on;
    plot(maxRefl,'g');
    plot(reflClean,'m');
    plot(numClean,'r');
    ylim([0 100]);
    hold off;
    drawnow;
    
    reflClean(badIndNum)=nan;
    meanRefl=mean(reflClean,'omitnan');
    
    %Radar times in UTC
    hh = floor(angTime/3600);
    mm = floor(mod(angTime,3600)/60);
    ss = floor(mod(angTime,60));  % this creates fractional seconds
    ff = floor(mod(angTime,1)*100);
    
    %
    % for i=1:length(badIndRefl);
    % output=['Bad reflectivity value of ',num2str(outlRefl(i)),' dBZ detected at ',...
    %     num2str(hh(badIndRefl(i))),':',num2str(mm(badIndRefl(i))),':',num2str(ss(badIndRefl(i)))];
    % disp(output);
    % end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate HCR bias
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    soundFile=strcat('/h/eol/romatsch/data/radar/soundings/D',radarFile(7:21),'_P.QC.eol');
    
    [attLiebe,attITU, attITU0, attITUW]=f_atten_layers(soundFile,f);
    
    c=3*10^8;
    tau=2.56*10^(-7);
    Ksquare=0.69;
    lambda=c/(f*10^9);
    
    dBsigma0=5.85;
    
    frac=(c*pi^5*tau*Ksquare)/(10^18*2*lambda^4);
    
    dBbiasITU=dBsigma0-meanRefl-2*attITU+10*log10(cos(deg2rad(radarAngle)))-10*log10(frac);
    dBbiasLiebe=dBsigma0-meanRefl-2*attLiebe+10*log10(cos(deg2rad(radarAngle)))-10*log10(frac);
    
    attLiebeAng=attLiebe/cosd(radarAngle);
    attITUAng=attITU/cosd(radarAngle);
    
    disp(strcat({'Mean reflectivity (dBZ): '},num2str(meanRefl)));
    disp('One way zenith attenuation:');
    disp(strcat({'Liebe: '},num2str(attLiebe)));
    disp(strcat({'ITU-R: '},num2str(attITU)));
    disp('HCR bias in dB:');
    disp(strcat({'Liebe: '},num2str(dBbiasLiebe)));
    disp(strcat({'ITU-R: '},num2str(dBbiasITU)));
    
    output(j,1)=meanRefl;
    output(j,2)=attLiebe;
    output(j,3)=attITU;
    output(j,4)=meanRefl+2*attLiebeAng;
    output(j,5)=meanRefl+2*attITUAng;
    output(j,6)=dBbiasLiebe;
    output(j,7)=dBbiasITU;
    output(j,8)=sum(~isnan(reflClean));
    
end

outTable = array2table(round(output,2), 'RowNames', allFileList);
outTable.Properties.VariableNames={'Refl' 'Att0Liebe'...
    'Att0ITU' 'TotReflLiebe' 'TotReflITU'...
    'dBbiasLiebe' 'dBbiasITU' 'NumData'};

%writetable(outTable,outFile,'WriteRowNames',true, 'Delimiter','\t')

spread=max(output,[],1)-min(output,[],1);