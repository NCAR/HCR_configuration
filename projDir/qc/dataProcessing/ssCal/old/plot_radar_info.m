% Plots flight and HCR info to be used in determination of surface
% reflectivity.
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

clear all
close all

addpath('/h/eol/romatsch/codes/matlab/hcrCalib/functions/');
addpath('/h/eol/romatsch/codes/matlab/hcrCalib/subCodes/');

direc = '/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/20180126';
%direc = '/scr/eldora1/rsfdata/cset/hcr/qc/cfradial/velcorr/10hz';
%direc = '/scr/eldora1/rsfdata/cset/hcr/cfradial/moments/10hz';
%direc = '/scr/eldora1/rsfdata/cset/hcr/cfradial/moments/100hz';
set(0,'DefaultTextInterpreter','none');

f_name = f_select_files(direc, '/cfrad*SUR.nc');

% make the assumption that Attributes will not change over the next
% few volumes (this is a serious assumption!)

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

klim = size(f_name,1);
fprintf(' %d files have been selected.\n',klim);

for kk=1:klim;
   fprintf('On file %d (of %d)\r',kk,klim);
   SWP = f_get_and_scale_hdf5_data(deblank(f_name(kk,:)), ATT,[ ...
       fields(ATT)' ]);
   if( kk == 1);
       dstr = [ char(SWP.time_coverage_start(:))' ];
       dstr = [ regexprep(dstr,'[A-Z]',' ') ];
       dstr = dstr(1:10);
   end;
   ngates = size(SWP.DBZ,1);
   nbeams = size(SWP.DBZ,2);

   PLT.rota  = [ PLT.rota;  SWP.rotation ];
   PLT.alt   = [ PLT.alt;   SWP.altitude ];
   PLT.pitch = [ PLT.pitch; SWP.pitch ];

   clear tmpsec;
   tmpsec = f_get_cfradnc_daysecs( SWP );
   
   startTimeIn=ncread(f_name(kk,:),'time_coverage_start')';
    startTime=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeRead=ncread(f_name(kk,:),'time')';
    
    PLT.time = [ PLT.time startTime+seconds(timeRead)];
   
   %PLT.time = [ PLT.time;  tmpsec ];
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
   
end;
fprintf('\nDone with files %d\n',kk);

plot_data

