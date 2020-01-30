function [PLT dB_bias_lb dB_bias_ITU avg_hdg avg_alt bias_lbITU_vec] = f_determine_bias(uniqueCases,filedir,indir,attITU,attLiebe,frq,bias_lbITU_vec)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
c = 3.0e8;
tau = 2.56e-7;
lambda = c/frq;
w_dielec_sq = .69;

temp = c * pi^5 * tau * w_dielec_sq /( 2*lambda^4 * 1e18);
rc_terms = 10 * log10(temp);

theta = 10.0;  % optimum angle for sigma0 determination
dB_sigma10 = 6.0;   % expected sigma0 at 10 degrees

clear PLT;
PLT = struct();
PLT.time = [];
PLT.rota = [];
PLT.refl = [];  % max reflectivity
PLT.vpwr = [];
PLT.alt  = [];
PLT.vel  = [];
PLT.pitch = [];
PLT.elev = [];  % elevation angle is referenced to horiz plane;
% pitch and roll corrected.  -80 = 170 or 190 rotation
PLT.angdif = []; % abs(rotation - 180) - (90 + elevation) (for downward)
PLT.roll = [];
PLT.range = [];
PLT.lon = [];
PLT.lat = [];
PLT.hdg = [];
PLT.sig0 = [];
PLT.sig0noBias = [];

oneCaseFileList=importdata([filedir,uniqueCases.casefiles{1}]);
for kk=1:size(oneCaseFileList,1);  % step through each input file in the given set
        
        indataOrig=[indir,oneCaseFileList{kk}(1:30)];
        
        getInfile=dir([indataOrig '*']);
        
        indata=[getInfile.folder,'/',getInfile.name];
        
        startTimeIn=ncread(indata,'time_coverage_start')';
        startTime=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(indata,'time')';
        
        timeIn=startTime+seconds(timeRead);
        
        rotaIn=ncread(indata,'rotation');
        pitchIn=ncread(indata,'pitch');
        rollIn=ncread(indata,'roll');
        elevIn=ncread(indata,'elevation');
        rangeIn=ncread(indata,'range');
        velIn=ncread(indata,'VEL');
        reflIn=ncread(indata,'DBZ');
        pwrIn=ncread(indata,'DBMVC');
        
        PLT.alt   = [ PLT.alt;   ncread(indata,'altitude')];
        PLT.hdg   = [ PLT.hdg;   ncread(indata,'heading') ];  % hope for no heading near North
        PLT.lat   = [ PLT.lat;   ncread(indata,'latitude') ];
        PLT.lon   = [ PLT.lon;   ncread(indata,'longitude') ];
        
        PLT.time = [ PLT.time;  timeIn' ];
        
        % There are high reflectivities in the first gates that need to be removed.
        badRow=8;
        empty=1;
        
        while (badRow<30 && empty>0)
            empty=length(find(~isnan(reflIn(badRow,:))));
            badRow=badRow+1;
        end;
        
        %sort out scans with too many reflectivity values
        dataSize=sum(~isnan(reflIn),1);
        tooMuchInd=find(dataSize>50);
        
        if badRow==30
            disp(['First row with only nans = ',num2str(badRow)]);
        else        
            reflIn(1:badRow-1,:)=nan;
        end
        
        [maxRefl maxGate]=nanmax(reflIn,[],1);
        maxRefl(tooMuchInd)=nan;
        
        %sort out scans with reflectivity values that are too high or low
        outRangeInd=find(maxRefl<0 | maxRefl>55);
        maxRefl(outRangeInd)=nan;
        
        %sort out scans where roll angle is too large or small
        outRollInd=find(rollIn<-2 | rollIn>2);
        maxRefl(outRollInd)=nan;
        
        maxGateLin=sub2ind(size(reflIn),maxGate,1:size(reflIn,2));
        emptyMax=find(isnan(maxRefl));
        
%         if max(rotaIn)>215 | min(rotaIn)<145
%             disp('Rotation out of range.');
%         end
        
        elev = atand(sqrt((tand(pitchIn)).^2 + ...
            (tand(rotaIn+ rollIn)./cosd(pitchIn)).^2));
        
        eldif = elevIn - elev + 90;
        
        maxpwr = pwrIn(maxGateLin);
        rangeMat=repmat(rangeIn,1,size(reflIn,2));
        sfcrng = rangeMat(maxGateLin);
        sfcvel = velIn(maxGateLin);
        
        if ~isempty(emptyMax)
            maxpwr(emptyMax)=nan;
            sfcrng(emptyMax)=nan;
            sfcvel(emptyMax)=nan;
            elev(emptyMax)=nan;
            eldif(emptyMax)=nan;
        end
        
        PLT.refl = [ PLT.refl; maxRefl' ];
        PLT.vel  = [ PLT.vel;  sfcvel'];
        PLT.vpwr = [ PLT.vpwr; maxpwr' ];
        PLT.range = [ PLT.range; sfcrng' ];
        PLT.elev  = [ PLT.elev; elev ];
        PLT.angdif = [ PLT.angdif; eldif ];
        PLT.rota = [ PLT.rota; rotaIn ];
        PLT.pitch = [ PLT.pitch; pitchIn ];
        PLT.roll = [ PLT.roll; rollIn ];
        
    end;
    
    timeInds=[];
    for ll=1:size(uniqueCases,1)
        % take only data in right time span
        begTime=datetime(str2num(uniqueCases.timest{ll}(1:4)),str2num(uniqueCases.timest{ll}(5:6)),str2num(uniqueCases.timest{ll}(7:8)),...
            str2num(uniqueCases.timest{ll}(10:11)),str2num(uniqueCases.timest{ll}(12:13)),str2num(uniqueCases.timest{ll}(14:15)));
        endTime=datetime(str2num(uniqueCases.timend{ll}(1:4)),str2num(uniqueCases.timend{ll}(5:6)),str2num(uniqueCases.timend{ll}(7:8)),...
            str2num(uniqueCases.timend{ll}(10:11)),str2num(uniqueCases.timend{ll}(12:13)),str2num(uniqueCases.timend{ll}(14:15)));
        
        timeInds=[timeInds; find(PLT.time>begTime+1/24/60/60 & PLT.time<endTime)];
    end
    
    fieldsPLT=fields(PLT);
    
    for jj=1:length(fieldsPLT)
        if ~isempty(PLT.(fieldsPLT{jj}))
            PLT.(fieldsPLT{jj})=PLT.(fieldsPLT{jj})(timeInds);
        end
    end
    %fprintf('\nDone with files %d\n',kk);
    clear A;
    % find beams near 10-deg elevation angle:
    A = find( abs(theta - PLT.elev) < 0.3);
    %average in non dB space
    [avg_dz,~,~]=dB_meanStd(PLT.refl(A));
    avg_alt = nanmean(PLT.alt);
    avg_hdg = nanmean(PLT.hdg);
    
    if ~isnan(attLiebe)
        dB_bias_lb = dB_sigma10 - 2*attLiebe  +  ...
            10*log10(cosd(PLT.elev(A))) - rc_terms - avg_dz;
        
        dB_bias_ITU = dB_sigma10 - 2*attITU +  ...
            10*log10(cosd(PLT.elev(A))) - rc_terms - avg_dz;
        
        %bias around 10 deg elevation angle for each single scan
        dB_bias_lb_vec= dB_sigma10 - 2*attLiebe  +  ...
            10*log10(cosd(PLT.elev(A))) - rc_terms - PLT.refl(A);
        
        dB_bias_ITU_vec= dB_sigma10 - 2*attITU  +  ...
            10*log10(cosd(PLT.elev(A))) - rc_terms - PLT.refl(A);
        
        %average bias around 10 deg for whole scanning event
%         [dBL,dBstdUL,dBstdDL]=dB_meanStd(dB_bias_lb_vec);
%         [dBITU,dBstdUITU,dBstdDITU]=dB_meanStd(dB_bias_ITU_vec);
        
        bias_lbITU_vec=cat(1,bias_lbITU_vec,cat(2,dB_bias_lb_vec,dB_bias_ITU_vec));
        
        %calculate sig0 measured without bias for all data, not just around
        %10 deg
        PLT.sig0noBias=PLT.refl + 2*attLiebe - 10*log10(cosd(PLT.elev)) + rc_terms;
        
        %PLT.sig0 = PLT.sig0noBias + dB_bias_lb;
        
    else
        PLT.sig0=nan;
        dB_bias_lb=nan;
        dB_bias_ITU=nan;
        PLT.sig0noBias=nan;
    end
    
end

