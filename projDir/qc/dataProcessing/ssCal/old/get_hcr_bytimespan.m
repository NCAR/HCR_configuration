function f_name = get_hcr_bytimespan( start_time, end_time, subdir1, subdir2, subdirs3);

% Uses a start_time and end_time to find the full list of path/file
% names for all scan types under the directory structure used for
% CfRadial data.  This routine uses get_sorted_file_list to find
% a full day of scans, concatenates days, then trims the start and end
% of the file list to match the desired start/end times.
%
% The routine requires that the start/end times are strings in the
% form 'ccyymmdd_hhmmss'.  Multiple subdirectory levels are
% specified, a holdover from multi-source comparison work.  The
% last subdir is a cell vector, allowing combination of sources
% such as 'sur', 'ppi', 'rhi', common last subdirs for S-Pol data.
%
% Remember that, within the RSF scheme, scans are often grouped first 
% by type of data (moments, covars, or whatever), then by 
% s or x band, then by scan type.
%
%  f_name = get_hcr_bytimespan('20110411_120000','20110411_175900',
%  'moments', 'sband');
%
% or 
%  f_name = get_hcr_bytimespan('20110411_120000','20110411_175900',
%  'moments', 'sband', {'sur'; 'sec'});
%
% Note that the get_sorted_file_list routine requires that a global
% root directory be set, but that this routine does not require/use
% that global variable.

% global rt_direc;   rt_direc = '/export/d2/rilling/DYN_test/cfradial';

% This routine requires our more modern cfradial file name convention, 
% where the file name starts with cfrad, and includes both the
% start and end time of the data within the file name. 

% Sanity checks on input times:

%start_time = char(start_time)
%end_time = char(end_time)

if( length(start_time) ~= 15 || length(end_time) ~= 15 );
    fprintf(['Input Time values do not have correct length in ' ...
            'get_file_list_bytimespan\n']);
    fprintf('Aborting.\n');
    return;
end;

entm = datenum(end_time,'yyyymmdd_HHMMSS');
sttm = datenum(start_time,'yyyymmdd_HHMMSS');

if( sttm >= entm );
    fprintf('Start time is greater than or equal to end time.\n');
    fprintf('Aborting\n');
    return;
end;

start_date = start_time(1:8);
end_date   = end_time(1:8);

f_name = [];
t_name = [];
my_date = start_date;
%class(my_date)

while ( my_date <= end_date );
   fn = get_sorted_file_list( char(my_date), subdir1, subdir2, subdirs3 );
   t_name  = [ t_name; fn ];
   my_date = datestr((datenum(my_date,'yyyymmdd') + 1), 'yyyymmdd');
end;
% assumes cfrad files with dual date/time (start and end) naming convention
mypos = regexp(t_name,'[1,2]\d\d\d\d\d\d\d_\d\d\d\d'); % find position of all date  strings in file names

s_time = [];
e_time = [];

for ii=1:numel(t_name);
    s_time = [ s_time datenum(fn{ii}(mypos{ii}(1):mypos{ii}(1)+ ...
                                     14),'yyyymmdd_HHMMSS') ];
    e_time = [ e_time datenum(fn{ii}(mypos{ii}(2):mypos{ii}(2)+ ...
                                     14),'yyyymmdd_HHMMSS')];
end;

A = find( s_time <= sttm, 1, 'last' );
B = find( e_time >= entm, 1, 'first' );

f_name = t_name(A:B);


