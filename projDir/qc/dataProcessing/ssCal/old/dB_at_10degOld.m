%Find the ocean surface reflectivity at 10 degree beam angle


clear all;
close all;

addpath('/h/eol/romatsch/matlab/radar/functions/');
addpath('/h/eol/romatsch/matlab/radar/subCodes/');

%Rotation angle (degree) at which reflectivity is wanted
degDB=10;
degPlusMinus=0.3; %Tolerance of readar angle

%Time tolerance in seconds
timePlusMinus=300;

%...read file list
radarFile='sound_20150722_193012.dat';
%radarFile='sound_20150724_165608.dat';
%radarFile='sound_20150724_230006.dat';
%radarFile='sound_20150801_214509.dat';
%radarFile='sound_20150812_201151.dat';

%radarFile='sound_20150709_222840.dat';
%radarFile='sound_20150709_224819.dat';
%radarFile='sound_20150803_170007.dat';
%radarFile='sound_20150807_193958.dat';
%radarFile='sound_20150807_212000.dat';

radarFileDir='/h/eol/romatsch/matlab/radar/fileLists/radarScans/';
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

totRota=[];
totTime=[];
totRefl=[];

for i=1:length(fileList);
    f_name = strcat(direc,soundDate,'/',cell2mat(fileList(i)));
    
    [ ATT GATT ] = f_get_hdf5_param_atts( f_name(1,:) );
    Flds = f_determine_radar_fields( ATT );
    % sFlds = select_params(Flds);
    % tFlds = determine_time_depends( SWP );
    
    if( isfield(GATT,'instrument_name') && size(GATT.instrument_name,2) > 0 )
        plat_name = GATT.instrument_name;
    else
        plat_name = 'UNKNOWN';
    end;
    
    clear PLT;
    PLT = struct();
    PLT.time = [];
    PLT.rota = [];
    PLT.refl = [];  % max reflectivity
    PLT.alt  = [];
    PLT.pitch = [];
    PLT.range = [];
    
    
    SWP = f_get_and_scale_hdf5_data(deblank(f_name), ATT,[ ...
        fields(ATT)' ]);
    dstr = [ char(SWP.time_coverage_start(:))' ];
    dstr = [ regexprep(dstr,'[A-Z]',' ') ];
    dstr = dstr(1:10);
    
    ngates = size(SWP.DBZ,1);
    nbeams = size(SWP.DBZ,2);
    
    PLT.rota  = [ PLT.rota;  SWP.rotation ];
    PLT.alt   = [ PLT.alt;   SWP.altitude ];
    PLT.pitch = [ PLT.pitch; SWP.pitch ];
    
    clear tmpsec;
    tmpsec = f_get_cfradnc_daysecs( SWP );
    PLT.time = [ PLT.time;  tmpsec ];
    if(size(tmpsec,1) ~= nbeams);
        fprintf('error on file: nbeams not equal to number of time points\n');
    end;
    
    SWP.DBZ(1:16,:) = NaN;  % there are high reflectivities in the
    % first gates that need to be removed.
    for jj=1:nbeams;
        clear tmpdz;
        if( SWP.rotation(jj) < 215.0 && SWP.rotation(jj) > 145.0)
            % smoothing is resource intensive.  try without.
            %           tmpdz = smooth( SWP.DBZ(:,jj),3,'lowess');
            %           A = find(isnan(SWP.DBZ(:,jj)));
            %           tmpdz(A) = NaN;
            tmpdz = SWP.DBZ(:,jj);
            maxgate = find( tmpdz == nanmax(tmpdz),1);
            if( ~(isnan(maxgate)) )
                maxdz = tmpdz(maxgate);
                sfcrng = SWP.range(maxgate);
            else
                maxdz = NaN;
                sfcrng = NaN;
            end;
        else;
            maxdz = NaN;
            sfcrng = NaN;
        end;
        PLT.refl = [ PLT.refl; maxdz ];
        PLT.range = [ PLT.range; sfcrng ];
    end;
    totRota=cat(1,totRota,PLT.rota);
    totTime=cat(1,totTime,PLT.time);
    totRefl=cat(1,totRefl,PLT.refl);
end;

plus10=find(totRota>(180+degDB-degPlusMinus) & totRota<(180+degDB+degPlusMinus));
minus10=find(totRota>(180-degDB-degPlusMinus) & totRota<(180-degDB+degPlusMinus));
goodInd=cat(1,plus10,minus10);
goodAngles=totRota(goodInd);

soundTimeLong=str2num(soundHH)*3600+str2num(soundMM)*60+str2num(soundSS);

timeAng=totTime(goodInd);
diffTime=abs(timeAng-soundTimeLong);
bestInd=find(diffTime<=timePlusMinus);
bestTime=timeAng(bestInd);
%veryBestInd=find(diffTime==min(diffTime));

hh = floor(bestTime/3600);
mm = floor(mod(bestTime,3600)/60);
ss = floor(mod(bestTime,60));  % this creates fractional seconds
ff = floor(mod(bestTime,1)*100);

goodRefl=totRefl(goodInd);
bestRefl=goodRefl(bestInd);
[reflClean, badIndRefl,outlRefl]=f_deleteoutliers(bestRefl,0.01);
meanRefl=mean(reflClean);
%veryBestRefl=goodRefl(veryBestInd)

%outAngle=abs(goodAngles(veryBestInd)-180)

output1=['Mean reflectivity ',num2str(meanRefl),' dBZ.'];
disp(output1);
for i=1:length(badIndRefl);
output2=['Bad reflectivity value of ',num2str(outlRefl(i)),' dBZ detected at ',...
    num2str(hh(badIndRefl(i))),':',num2str(mm(badIndRefl(i))),':',num2str(ss(badIndRefl(i)))];
disp(output2);
end;
