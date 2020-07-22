% Plots flight and HCR info to be used in determination of surface
% reflectivity.  Determines bias according to Rilling technique.
% Uses a calibration input file, with attenuations pre-computed
% by Romatschke.
%
% Creates multi-panel plots over a selected, contiguous, number of
% data files.  Panels include:
%
%   -- scan rotation angle (rotation is in a plane perpendicular
%                           to flight direction)
%   -- max reflectivity (when pointing down, only)
%   -- range to max reflectivity (when pointing down)
%   -- Altitude
%   -- Pitch
%
% all parameters are plotted vs time, where time is in seconds from
% the start of the first file, or seconds from the start of the day.
%
% Note that for sea surface scans, we operate in vertical receive, only.

clear all;
close all;

addpath('/h/eol/romatsch/matlab/radar/functions/');

figdir='/h/eol/romatsch/radar/figures/cases/';

indir='/scr/eldora1/rsfdata/aristo-17/hcr/cfradial/moments/10hz/20170302/';

fileList={{'cfrad.20170302_214800.132_to_20170302_214900.046_HCR_SUR.nc';...
    'cfrad.20170302_214900.147_to_20170302_215000.060_HCR_SUR.nc';...
    'cfrad.20170302_215000.162_to_20170302_215100.075_HCR_SUR.nc';...
    'cfrad.20170302_215100.176_to_20170302_215200.089_HCR_SUR.nc';...
    'cfrad.20170302_215200.191_to_20170302_215300.003_HCR_SUR.nc';...
    'cfrad.20170302_215300.104_to_20170302_215400.017_HCR_SUR.nc';...
    'cfrad.20170302_215400.119_to_20170302_215500.032_HCR_SUR.nc';...
    'cfrad.20170302_215500.133_to_20170302_215600.046_HCR_SUR.nc';...
    'cfrad.20170302_215600.148_to_20170302_215700.061_HCR_SUR.nc'}};

fileList{end+1}={'cfrad.20170302_203600.197_to_20170302_203700.009_HCR_SUR.nc';...
    'cfrad.20170302_203700.110_to_20170302_203800.023_HCR_SUR.nc';...
    'cfrad.20170302_203800.125_to_20170302_203900.038_HCR_SUR.nc';...
    'cfrad.20170302_203900.139_to_20170302_204000.052_HCR_SUR.nc';...
    'cfrad.20170302_204000.154_to_20170302_204100.067_HCR_SUR.nc';...
    'cfrad.20170302_204100.168_to_20170302_204200.082_HCR_SUR.nc';...
    'cfrad.20170302_204200.183_to_20170302_204300.096_HCR_SUR.nc';...
    'cfrad.20170302_204300.198_to_20170302_204400.010_HCR_SUR.nc';...
    'cfrad.20170302_204400.111_to_20170302_204500.024_HCR_SUR.nc';...
    'cfrad.20170302_204500.126_to_20170302_204600.039_HCR_SUR.nc'};

fileList{end+1}={'cfrad.20170302_213200.102_to_20170302_213300.015_HCR_SUR.nc';...
    'cfrad.20170302_213300.116_to_20170302_213400.030_HCR_SUR.nc';...
    'cfrad.20170302_213400.131_to_20170302_213500.044_HCR_SUR.nc';...
    'cfrad.20170302_213500.146_to_20170302_213600.059_HCR_SUR.nc';...
    'cfrad.20170302_213600.160_to_20170302_213700.073_HCR_SUR.nc';...
    'cfrad.20170302_213700.175_to_20170302_213800.088_HCR_SUR.nc';...
    'cfrad.20170302_213800.189_to_20170302_213900.001_HCR_SUR.nc';...
    'cfrad.20170302_213900.103_to_20170302_214000.016_HCR_SUR.nc';...
    'cfrad.20170302_214000.117_to_20170302_214100.030_HCR_SUR.nc'};

fileList{end+1}={'cfrad.20170302_220200.134_to_20170302_220300.047_HCR_SUR.nc';...
    'cfrad.20170302_220300.149_to_20170302_220400.062_HCR_SUR.nc';...
    'cfrad.20170302_220400.163_to_20170302_220500.076_HCR_SUR.nc';...
    'cfrad.20170302_220500.178_to_20170302_220600.091_HCR_SUR.nc';...
    'cfrad.20170302_220600.192_to_20170302_220700.004_HCR_SUR.nc';...
    'cfrad.20170302_220700.106_to_20170302_220800.019_HCR_SUR.nc';...
    'cfrad.20170302_220800.120_to_20170302_220900.033_HCR_SUR.nc';...
    'cfrad.20170302_220900.135_to_20170302_221000.048_HCR_SUR.nc';...
    'cfrad.20170302_221000.149_to_20170302_221100.063_HCR_SUR.nc';...
    'cfrad.20170302_221100.164_to_20170302_221200.077_HCR_SUR.nc';...
    'cfrad.20170302_221200.179_to_20170302_221300.092_HCR_SUR.nc'};

fileList{end+1}={'cfrad.20170302_222800.108_to_20170302_222900.021_HCR_SUR.nc';...
    'cfrad.20170302_222900.122_to_20170302_223000.036_HCR_SUR.nc';...
    'cfrad.20170302_223000.137_to_20170302_223100.050_HCR_SUR.nc';...
    'cfrad.20170302_223100.152_to_20170302_223200.065_HCR_SUR.nc';...
    'cfrad.20170302_223200.166_to_20170302_223300.079_HCR_SUR.nc';...
    'cfrad.20170302_223300.181_to_20170302_223400.094_HCR_SUR.nc';...
    'cfrad.20170302_223400.195_to_20170302_223500.007_HCR_SUR.nc'};

global h;  % used in axis scaling/resizing
global rt_direc;
% non-qc directories:
% direc = '/scr/eldora1/rsfdata/cset/hcr/cfradial/moments/10hz';
% rt_direc = '/scr/eldora1/rsfdata/cset/hcr/';

% qc directories
direc = '/scr/eldora1/rsfdata/cset/hcr/qc/cfradial/moments/10hz';
rt_direc = '/scr/eldora1/rsfdata/cset/hcr/qc/';

set(0,'DefaultTextInterpreter','none');

% recover the calibration info from a text file
% (you need to know which fields are text and which are numeric)
% A = importdata('cal.temp',' ',1);
A = importdata('cal_ARISTO.text',' ',1);
np = numel(A.textdata(:,1));
nd = numel(A.data(:,1));
cal.(A.textdata{1,1}) = A.textdata(2:np,1);
cal.(A.textdata{1,2}) = A.textdata(2:np,2);
cal.(A.textdata{1,3}) = A.textdata(2:np,3);
cal.(A.textdata{1,4}) = A.textdata(2:np,4);
cal.(A.textdata{1,5}) = A.data(:,1);
cal.(A.textdata{1,6}) = A.data(:,2);

% set up radar parameters, and partial radar constant term

c = 3.0e8;
cf = pi/180;
frq = 9.4406e10;
tau = 2.56e-7;
lambda = c/frq;
w_dielec_sq = .69;

temp = c * pi^5 * tau * w_dielec_sq /( 2*lambda^4 * 1e18);
rc_terms = 10 * log10(temp);

theta = 10.0;  % optimum angle for sigma0 determination
dB_sigma10 = 6.0;   % expected sigma0 at 10 degrees

for ncal=1:nd;
    close all;
    
    fileListInd=fileList{ncal};
    
    f_name={};
    for ii=1:length(fileListInd)
        f_name{end+1} = [indir fileListInd{ii}];
    end
    
    f_name=f_name';
    
    [ ATT GATT ] = get_hdf5_param_atts( char(f_name(1,:)) );
    Flds = determine_radar_fields( ATT );
    
    if( isfield(GATT,'instrument_name') && size(GATT.instrument_name,2) > 0 )
        plat_name = GATT.instrument_name;
    else
        plat_name = 'UNKNOWN';
    end;
    
    % initialize accumlation arrays for each pass through main loop
    
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
    
    % system parameters appear once for each input cfrad file
    clear SYS;
    SYS = struct();
    SYS.EikTemp = [];
    SYS.CathVolt = [];
    SYS.XMTemp = [];
    SYS.RfDetPwr = [];
    SYS.PloTemp = [];
    SYS.VLnaTemp = [];
    SYS.RfDetTemp = [];
    %SYS.CorrVc = [];
    SYS.PsVolt = [];
    
    klim = size(f_name,1);
    for kk=1:klim;  % step through each input file in the given set
        %   fprintf('On file %d (of %d)\r',kk,klim);
        SWP = get_and_scale_hdf5_data(char(f_name(kk,:)), ATT,[ ...
            fields(ATT)' ]);
        if( kk == 1);
            dstr = [ char(SWP.time_coverage_start(:))' ];
            fdstr = [ regexprep(dstr,'[A-Z]',' ') ];
            dstr = fdstr(1:10);
        end;
        
        % tedious and repetitive: extract system info stored in SWP.status_xml
        xml = char([SWP.status_xml{:}]);
        
        % extract the Klystron temperature value from the XML info
        temp_str = extractBetween(xml,'<EikTemp>','</EikTemp>');
        temp_str=temp_str{:};
        if( length(temp_str) == 0 )
            temp = NaN;
        else
            temp = sscanf(temp_str, '%f');
        end;
        SYS.EikTemp = [ SYS.EikTemp; temp ];
        
        temp_str = extractBetween(xml,'<CathodeVoltage>','</CathodeVoltage>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.CathVolt = [ SYS.CathVolt; temp ];
        
        temp_str = extractBetween(xml,'<XmitterTemp>','</XmitterTemp>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.XMTemp = [ SYS.XMTemp; temp ];
        
        temp_str = extractBetween(xml,'<DetectedRfPower>','</DetectedRfPower>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.RfDetPwr = [ SYS.RfDetPwr; temp ];
        
        temp_str = extractBetween(xml,'<PloTemp>','</PloTemp>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.PloTemp = [ SYS.PloTemp; temp ];
        
        temp_str = extractBetween(xml,'<VLnaTemp>','</VLnaTemp>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.VLnaTemp = [ SYS.VLnaTemp; temp ];
        
        temp_str = extractBetween(xml,'<RfDetectorTemp>','</RfDetectorTemp>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.RfDetTemp = [ SYS.RfDetTemp; temp ];
        
        %     temp_str = extractBetween(xml,'<gainCorrectionVc>','</gainCorrectionVc>');
        %     temp_str=temp_str{:};
        %     if( length(temp_str) == 0 )
        %         temp = NaN;
        %     else
        %         temp = sscanf(temp_str, '%f');
        %     end;
        %     SYS.CorrVc = [ SYS.CorrVc; temp ];
        
        temp_str = extractBetween(xml,'<PsVoltage>','</PsVoltage>');
        temp_str=temp_str{:};
        temp = sscanf(temp_str, '%f');
        SYS.PsVolt = [ SYS.PsVolt; temp ];
        
        ngates = size(SWP.DBZ,1);
        nbeams = size(SWP.DBZ,2);
        % accumulate multi-sweep attitude information
        PLT.rota  = [ PLT.rota;  SWP.rotation ];
        PLT.alt   = [ PLT.alt;   SWP.altitude ];
        PLT.pitch = [ PLT.pitch; SWP.pitch ];
        PLT.roll  = [ PLT.roll;  SWP.roll ];
        PLT.hdg   = [ PLT.hdg;   SWP.heading ];  % hope for no heading near North
        PLT.lat   = [ PLT.lat;   SWP.latitude ];
        PLT.lon   = [ PLT.lon;   SWP.longitude ];
        
        clear tmpsec;
        tmpsec = get_cfradnc_daysecs( SWP );
        PLT.time = [ PLT.time;  tmpsec ];
        if(size(tmpsec,1) ~= nbeams);
            fprintf('error on file: nbeams not equal to number of time points\n');
        end;
        
        SWP.DBZ(1:16,:) = NaN;  % there are high reflectivities in the
        % first gates that need to be removed.
        for jj=1:nbeams;
            clear tmpdz; clear tmppwr;
            if( SWP.rotation(jj) < 215.0 && SWP.rotation(jj) > 145.0)
                % smoothing is resource intensive.  try without.
                %           tmpdz = smooth( SWP.DBZ(:,jj),3,'lowess');
                %           A = find(isnan(SWP.DBZ(:,jj)));
                %           tmpdz(A) = NaN;
                tmpdz = SWP.DBZ(:,jj);
                tmppwr = SWP.DBMVC(:,jj);
                maxgate = find( tmpdz == nanmax(tmpdz),1);
                %      elev = SWP.elevation(jj);
                elev = atand(sqrt( (tand(SWP.pitch(jj)))^2 + ...
                    (tand(SWP.rotation(jj)+ SWP.roll(jj))/cosd(SWP.pitch(jj)))^2));
                %      eldif = abs( SWP.rotation(jj) + SWP.roll(jj) - 180) - (90 + elev);
                eldif = SWP.elevation(jj) - elev + 90;
                if( ~(isnan(maxgate)) )
                    maxdz = tmpdz(maxgate);
                    maxpwr = tmppwr(maxgate);
                    sfcrng = SWP.range(maxgate);
                    sfcvel = SWP.VEL(maxgate,jj);
                else
                    maxdz = NaN;
                    maxpwr = NaN;
                    sfcrng = NaN;
                    sfcvel = NaN;
                end;
            else;
                maxdz  = NaN;
                maxpwr = NaN;
                sfcrng = NaN;
                sfcvel = NaN;
                elev   = NaN;
                eldif  = NaN;
            end;
            % accumulate information on a beam-by-beam basis
            PLT.refl = [ PLT.refl; maxdz ];
            PLT.vel  = [ PLT.vel;  sfcvel];
            PLT.vpwr = [ PLT.vpwr; maxpwr ];
            PLT.range = [ PLT.range; sfcrng ];
            PLT.elev  = [ PLT.elev; elev ];
            PLT.angdif = [ PLT.angdif; eldif ];
        end;
    end;
    %fprintf('\nDone with files %d\n',kk);
    clear A;
    % find beams near 10-deg elevation angle:
    A = find( abs(10 - PLT.elev) < 0.3);
    avg_refl_mW = nanmean(10.^(PLT.refl(A)/10.));
    avg_dz = 10*log10(avg_refl_mW);
    avg_alt = nanmean(PLT.alt);
    avg_hdg = nanmean(PLT.hdg);
    
    avg_EikTemp = nanmean(SYS.EikTemp);
    avg_CathVolt = nanmean(SYS.CathVolt);
    avg_XMTemp = nanmean(SYS.XMTemp);
    avg_RfDetPwr = nanmean(SYS.RfDetPwr); % value does not change
    avg_PloTemp = nanmean(SYS.PloTemp);
    avg_VLnaTemp = nanmean(SYS.VLnaTemp);
    avg_RfDetTemp = nanmean(SYS.RfDetTemp);
    %avg_CorrVc = nanmean(SYS.CorrVc);
    avg_PsVolt = nanmean(SYS.PsVolt);
    
    
    dB_bias_lb = dB_sigma10 - 2*cal.at_lieb(ncal)  +  ...
        10*log10(cosd(theta)) - rc_terms - avg_dz;
    
    dB_bias_ITU = dB_sigma10 - 2*cal.at_ITU(ncal) +  ...
        10*log10(cosd(theta)) - rc_terms - avg_dz;
    
    PLT.smth_dz = [];
    PLT.smth_dz = smooth(PLT.refl,5);
    % call the script to plot scan angle, reflectivity, aircraft
    % pitch/roll, etc.  (uses local variables names)
    plot_hcr_scan_data;
    plabel = [figdir 'ARISTO_CAL_' cal.timest{ncal} '.png'];
    saveas(f,plabel);
    
    % using Liebe bias, compute sigma0 for all elevation angles, and
    % plot (since Liebe attenuation is very near to ITU, did not do both)
    
    clear sisg0; sig0 = [];
    sig0 = PLT.refl + 2*cal.at_lieb(ncal) - 10*log10(cosd(PLT.elev)) + ...
        rc_terms + dB_bias_lb;
    
    %divide into rotation angles bigger or smaller than 180
    sig0smaller=sig0(find(PLT.rota<=180));
    elevSmaller=PLT.elev(find(PLT.rota<=180));
    
    sig0bigger=sig0(find(PLT.rota>180));
    elevBigger=PLT.elev(find(PLT.rota>180));
    
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    %plot(PLT.elev,sig0);
    hold on
    plot(elevSmaller,sig0smaller,'r');
    plot(elevBigger,sig0bigger,'b');
    set(gca,'XLim',[0,25]);
    set(gca,'YLim',[-20,20]);
    xlabel('Incidence Angle, degrees');
    ylabel('Normalized Radar Cross Section, dB');
    title(['Bias-corrected Sigma0 vs Incidence Angle ' fdstr]);
    set(0,'DefaultTextInterpreter','none');
    tstr = sprintf(['Sfc wind = %s \nAvg Hdg  = %5.0f deg\nAvg Alt = ' ...
        '%7.0f\nLiebe bias = %4.2f dB\n2way attn = %5.2f dB'], cal.sfcwind{ncal}, ...
        avg_hdg, avg_alt,dB_bias_lb, 2*cal.at_lieb(ncal));
    text( 15,12,tstr);
    legend('Rot angle <= 180 deg','Rot angle > 180 deg','location','southwest');
    plabel = [figdir 'ARISTO_SIGMA_' cal.timest{ncal} '.png'];
    saveas(f2,plabel);
    
end