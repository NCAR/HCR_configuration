#!/usr/bin/perl

use Getopt::Std;

#require "getopt.pl";
use Time::Local;

$ENV{TZ} = "GMT";
$odir = "";

use constant TRUE  => 1;     # preprocessor command to substitute 1 for TRUE in text
use constant FALSE => 0;

if( !getopts('d:f:F:I:m:N:P:s:S:t:T:aCiMV') ) {   # note that there is a problem passing a zero value as an opt

print <<_EOF_;

  USAGE: 
		$0 [OPTIONS] [FILES]
	    or
	        ls -1d *files* | $0 [OPTIONS]

	Assists in tar-ing up time-stamped data files, particularly 
	radar sweep files, whether in netCDF or DORADE format.  This
	routine may also be used for almost any time-stamped file, but
	results may vary, depending upon the actual meta information
	included in the filename.

	This routine builds a file that can then be executed to do the 
	actual tar work.  Review your tar file before executing!

	The input list of files, or the file listing itself, must be
	presented as a time-ordered list, with time increasing.

	Care is taken to ensure that a list of filename to tar can be
	input in multiple ways.  This can avoid problems with an input
	argument list being too long, and can allow careful filtering 
	of files for inclusion.

	Directory structures can be preserved, but again, actual results 
	may vary.

	This routine is primarily designed for radar data files, where the
	filename includes components for the timestamp, platform name,
	sweep type, volume number, and tilt angle.  In the case of other
	time-stamped files, default values are used for the radar-specific
	meta information.

  OPTIONS:

	-a   
	     tar all time-stamped files in this directory; in the case
	     of subdirectories, will also include time-stamped files
	     from the subdirectories

	-C
	     Show continuous periods of data files.  This option de-activates
	     any tar operations.  Allowed modifiers with this option are limited
	     to f, a, and t.  Scan type is ignored.  This option is included 
	     simply because all the I/O and processing required to tar the files, 
	     is the same as that required to find continuous periods.

	-d DIR
	     Declare an output directory, DIR

	-f filename
	     Declare an input file that contains a list of files to tar

	-F fmt
	     Short format descriptor (e.g., swp, nc, fof)

	-i
	     Ignore scan type

	-t nnn
	     Start a new tar archive after a time gap of nnn seconds
	     between files

	-I 'fname1 fname2 ...'
	     Space delimited, and may include wildcards; this
	     information will be transcribed exactly as the first part
	     of the tar statement.  Include the named non-sweep files
	     in each of the tar archives; for this option, use a full
	     relative pathname for the files.  This option is most
	     often used to include a README file in all archives.

  WAYS TO BUNDLE FILES

	-m nn
	     Create a new archive every nn minutes, where mod(60,nn) = 0,
	     as well as some other values.  Allowed values of nn are 
	     1,2,3,4,5,6,10,12,15,20,30,60, 120, 180, 240, 360, 480, 
	     720, and 1440.
	     
	     Selection of the m option causes other bundling options to 
	     be ignored.

	-M 
	     tar by multiple scan volumes, regardless of scan type,
	     until a maximum size or maximum period length is reached
	     (the small s and the T options); may not need this if 
	     the -i option is used.

	-N   explicitly declare the platform name

	-S nn
	     Exclusively tar by desired archive filesize, in MegaBytes
	     (ignore scan type, time gaps, or anything else)

	-s nn
	     tar by desired archive filesize, in MegaBytes, generally
	     in conjunction with other switches, such as scan type 
	     or time interval.  Preferentially use the -S option (the other 
	     option!) whenever you have a well-structured, continuous, 
	     non-radar data set.

	-T nn
	     tar by maximum time period length, in minutes.  This option
	     incorporates use of -t and volume information.

	-V
	     Use the scan volume number (obviously will not work for
	     files that are not named by scan file conventions).  This
	     over-rides the default of trying to determine a volume by
	     reviewing the fixed angle.  This option just changes how
	     the volume is defined, and does not alter requested
	     behavior based on archive size or changes in scan type.
	     Use this option when working with RHIs, or where the
	     elevation angle may dip to a low scan in the middle of a
	     scan volume.  (Alternatively, consider this as "ignore
	     fixed angle")

  NOTES
	It is unlikely that all options, or all combinations of options, 
	will work under all circumstances.  Be sure to inspect your job
	file before exectuting to ensure that the archive process behaves
	correctly.

	In the special case of ELDORA files, you may want to archive 
	TA and TF ("tail-aft" and "tail-forward", respectively) files
	separately.  If so, you should run this tar job twice, or start
	with the TA and TF files in separate directories.

  EXAMPLES:

	When working with files from multiple directories, the following
	call will fail:

		$0 -t 70 201003*

	This will work:

		$0 -t 70 201003*/PLATFORM_*

_EOF_
exit;
}

$use_size_only = FALSE;
$use_size = FALSE;
$use_time = FALSE;
$use_minutes = FALSE;
$show_continuous = FALSE;
$use_multiv = FALSE;
$use_vol_num = FALSE;

if ( $opt_d ) {         	# set an output directory
    $odir = $opt_d;
    if( substr($odir,length($odir)-1,1) ne "/" ) {
	$odir .= "/";
    }
}
if ( $opt_m ) {         	# opt_m does not exist if a zero was passed with -m 0
    $m_int = $opt_m;
    $tststr = sprintf("_%s_",$opt_m);

    if ( ($m_int > 0) && index("_1_2_3_4_5_6_10_12_15_20_30_60_120_180_240_360_480_720_1440_",$tststr) >= 0 ) {
	$sec_int = $m_int * 60;
	$use_minutes = TRUE;
    } else {
	printf("\nError.  m must divide 60 without a remainder, or be an allowed fraction of a day.\n");
	printf("Allowed values of m are 1,2,3,4,5,6,10,12,15,20,30,60,\n120, 180, 240, 360, 480, 720, and 1440.\n\n");
	exit;
    }
}

$fsize = 200;       # fsize default in MBytes; could be some problems with setting this, here.

if ( $opt_S ) {
    $fsize = $opt_S;    # set appoximate average tar archive size, fsize, in MBytes
    $use_size_only = TRUE;
}

if ( $opt_s ) {
    $fsize = $opt_s;    # set appoximate tar archive size, fsize, in MBytes
    $use_size = TRUE;
}


if ( $opt_F ) {		# explicitly include a format tag in the archive filename
    $fmt = $opt_F;
    $fmt = "_" . $fmt;
}
else {
    $fmt = "";
}

if ( $opt_I ) {
    $incl = $opt_I;
    printf("got an include set:  +++%s+++\n", $incl);
}
else {
    $incl = "";
}

if ( $opt_t ) {		# set an allowed time gap used to trigger a new tar archive
    $tgap = $opt_t;
}
else {
    $tgap = 100;
}


if ( $opt_T ) {		# set a time limit (in minutes) on the max length of an archive
    			# times start at non-uniform intervals (see -m option for contrast)
    $archtime = $opt_T;
    $use_time = TRUE;
}
else {
    $archtime = 600;	# default, but will not work because $use_time is FALSE
}

if ( $opt_N ) {
    $thisname = $opt_N;
} else {
    $thisname = "";
}


if ( $opt_P ) {
    $my_plat = TRUE;
    $this_plat = $opt_P;
}

if ( $opt_a ) {
    $all = TRUE;
}

if ( $opt_i ) {
    $ignore_type = TRUE;
} else {
    $ignore_type = FALSE;
}

if ( $opt_C ) {
    $show_continuous = TRUE;
}

if ( $opt_M ) {
    $use_multiv = TRUE;
}

if ( $opt_V ) {
    $use_vol_num = TRUE;
}

# four different ways to get a list of files to process:
#    as a list from an input file, as shown with the -f switch
#    as "all", with the -a switch
#    as the remaining ARGVs
#    as piped filenames
# the particular methods are mutually exclusive

if ( $opt_f ) {
    @fnames = `cat $opt_f`;
    for($i=0; $i<=$#fnames; $i++) { $fnames[$i] =~ chomp($fnames[$i]);}	
} elsif ( $all == TRUE ) {
    @fnames = `ls -1`;
    for($i=0; $i<=$#fnames; $i++) { $fnames[$i] =~ chomp($fnames[$i]);}	
} elsif( @ARGV ) {
    for(@ARGV) {
	push @fnames, $_;
    }
} else {
    while(<>) {
	chomp;
	push @fnames, $_;
    }
}

#print "$fnames[$#ARGV] \n";

# allows matching of dates/times in certain formats.
# Currently planned to work with files using the following conventions:
#
#       ncswp_NAME_ccyymmdd_hhmmss.*.nc (NAME can be any character string)
#	cfrad.ccyymmdd_hhmmss_*.nc  (no system name available)
#       NAME_ccyymmddhhmmss.*  or NAME.ccyymmddhhmmss.*
#	NAME_ccyymmdd_hhmmss.* or NAME.ccyymmdd_hhmmss.*
#       swp.1yymmddhhmmss.NAME*
#       swp.yymmddhhmmss.NAME*
#	(and a few other logical variations)
#
#Do not mix-and-match different conventions within the same call to this command.
#

push @sweeps, "swp.dummy";
push @tar_time, "";
push @swptime, 0;
push @s_size, 0;
push @fx_ang, 0;
push @vol_num, 0;
push @stype, "NUL";
push @name, "NUL";

ITER: for(@fnames) {
    chomp;
    $this_file = $_;
    $stemp = $this_file;
    if( ($stemp =~ m/\// ) ) { 		# look only at the filename component
	$stemp = substr($stemp,1+rindex($stemp,"/"));
    }

# Test for the most restrictive cases, first.  This will likely do a 
# better job of protecting against the unforeseen.  The most 
# restrictive are likely to be the dorade or netCDF sweepfiles, 
# followed by generic, dated data files.  Only the sweepfiles will 
# have all the various desired info.

# All allowed formats have deliminators (a . or _) between the
# date/time info and the rest of the filename components.  No
# exceptions.  All allowed formats require an 8-digit year (except
# swp) in ccyyddmm order, and a 6-digit time (time to nearest second).

# We could get fancier with the string matches, but showing explicit
# cases is not a bad thing, either.

    if ( ($stemp =~ m/swp\.1\d\d\d\d\d\d\d\d\d\d\d\d\./g) ) {	# DORADE sweepfiles, recent years
	$ind = pos $stemp;  # pos is the position just past the match
	$tststr = substr($stemp,$ind-13,12);
	($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/,$tststr,7));
	$yr += 2000;
        ($pref,$scan_t,$vanum) = split(/_/,$stemp,3);
        $vol = substr($vanum,1);
        $i = rindex($pref,".",length($pref)-3);
        $ang = substr($pref,$i+1);
	($junk,$dateinfo,$plat,$junk) = (split(/\./,$pref,4));
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d vol=%d ang=%f type=%s\n",$plat,$yr,$mon,$day,$hr,$min,$sec,$vol,$ang,$scan_t);
    } elsif ( ($stemp =~ m/swp\.\d\d\d\d\d\d\d\d\d\d\d\d\./g) ) {	# DORADE sweepfiles, pre 2000
	$ind = pos $stemp;
	$tststr = substr($stemp,$ind-13,12);
	($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/,$tststr,7));
	$yr += 1900;
        ($pref,$scan_t,$vanum) = split(/_/,$stemp,3);
        $vol = substr($vanum,1);
        $i = rindex($pref,".",length($pref)-3);
        $ang = substr($pref,$i+1);
	($junk,$dateinfo,$plat,$junk) = (split(/\./,$pref,4));
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d vol=%d ang=%f type=%s\n",$plat,$yr,$mon,$day,$hr,$min,$sec,$vol,$ang,$scan_t);
    } elsif ( ($stemp =~ m/ncswp/ ) && ($stemp =~ m/_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\./g ) ) {	# netCDF sweep files
	$ind = pos $stemp;
        $tststr = substr($stemp,$ind-16,15);
        ($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d\d\d)(\d\d)(\d\d)_(\d\d)(\d\d)(\d\d)/,$tststr,7));
        ($pref,$plat,$d1,$d2,$vanum,$scn,$ang,$scan_t,$suff) = split(/_/,$stemp,9);
        $vol = substr($vanum,1);
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d vol=%d ang=%f type=%s\n",$plat,$yr,$mon,$day,$hr,$min,$sec,$vol,$ang,$scan_t);
    } elsif ( ($stemp =~ m/cfrad/ ) && ($stemp =~ m/\.\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\_/g ) ) {	# CfRadial files
	$ind = pos $stemp;
        $tststr = substr($stemp,$ind-16,15);
        ($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d\d\d)(\d\d)(\d\d)_(\d\d)(\d\d)(\d\d)/,$tststr,7));
	($pref,$core,$suff,$nc) = split(/\./,$stemp,4);
	$plat = $thisname;  # for the current condition where cfrad files do not have a platform name
        ($d1,$d2,$msec,$vanum,$scn,$ang) = split(/_/,$core,6);
	$ang = substr($ang,2);  $ang = $ang . "." . substr($suff,0,2);
#	printf("%s\n%s\n%s\n%s\n",$msec,$vanum,$scn,$ang);
	($junk,$scan_t) = split(/_/,$suff,2);
        $vol = substr($vanum,1);
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d vol=%d ang=%s type=%s\n",$plat,$yr,$mon,$day,$hr,$min,$sec,$vol,$ang,$scan_t);
    } elsif ( ($stemp =~ m/_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\./g)   || 
	      ($stemp =~ m/_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d_/g)    ||
	      ($stemp =~ m/^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d_/g)    ||   # time series name format
	      ($stemp =~ m/^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\./g)   ||   # yet another time series name format
	      ($stemp =~ m/\.\d\d\d\d\d\d\d\d\.\d\d\d\d\d\d\./g) ||
	      ($stemp =~ m/\.\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\./g) ||
	      ($stemp =~ m/\.\d\d\d\d\d\d\d\dT\d\d\d\d\d\d[_.]/g) || 
	      ($stemp =~ m/\_\d\d\d\d\d\d\d\dT\d\d\d\d\d\d[_.]/g) )  {
	$ind = pos $stemp;
        $tststr = substr($stemp,$ind-16,15);
        ($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d\d\d)(\d\d)(\d\d)[._T](\d\d)(\d\d)(\d\d)/,$tststr,7));
	($plat,$junk) = (split(/[._]/,$stemp,2));
	$scan_t = "";
	$vol    = 0;
	$ang    = 0;
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d\n",$plat,$yr,$mon,$day,$hr,$min,$sec);
    } elsif ( ($stemp =~ m/_\d\d\d\d\d\d\d\d\d\d\d\d\d\d\./g) ||
	      ($stemp =~ m/_\d\d\d\d\d\d\d\d\d\d\d\d\d\d_/g)  ||
	      ($stemp =~ m/\.\d\d\d\d\d\d\d\d\d\d\d\d\d\d\./g)) {
	$ind = pos $stemp;
        $tststr = substr($stemp,$ind-15,14);
        ($junk,$yr,$mon,$day,$hr,$min,$sec) = (split(/(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/,$tststr,7));
	($plat,$junk) = (split(/[._]/,$stemp,2));
	$scan_t = "";
	$vol    = 0;
	$ang    = 0;
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d\n",$plat,$yr,$mon,$day,$hr,$min,$sec);
#   a somewhat abbreviated date format that does not include the seconds field:
    } elsif ( ($stemp =~ m/_\d\d\d\d\d\d\d\d\d\d\d\d\./g) ||
	      ($stemp =~ m/_\d\d\d\d\d\d\d\d\d\d\d\d_/g)  ||
	      ($stemp =~ m/\.\d\d\d\d\d\d\d\d\d\d\d\d\./g)) {
	$ind = pos $stemp;
        $tststr = substr($stemp,$ind-13,12);
#	printf("tststr = **%s**\n",$tststr);
        ($junk,$yr,$mon,$day,$hr,$min) = (split(/(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/,$tststr,6));
	$sec = 0;
	($plat,$junk) = (split(/[._]/,$stemp,2));
	$scan_t = "";
	$vol    = 0;
	$ang    = 0;
#        printf("plat=%s %d:%02d:%02d:%02d:%02d:%02d\n",$plat,$yr,$mon,$day,$hr,$min,$sec);
    } else {
        printf("INVALID date-based file name: %s\n", $stemp);
        next ITER;
    }
    $xtime = timegm($sec,$min,$hr,$day,$mon-1,$yr-1900);
    $temp = sprintf("%04d%02d%02d_%02d%02d%02d",$yr,$mon,$day,$hr,$min,$sec);
    push @tar_time, $temp;
    push @sweeps, $this_file;
    push @swptime, $xtime;
    push @s_size, (stat($this_file))[7]/1000000.0;
    push @fx_ang, $ang;
    push @vol_num, $vol;
    push @stype, $scan_t;
    push @name, $plat;

#    printf("unix_time = %d\n",$xtime);
}

$sweeps[0]  = $sweeps[1];    # for convenience, set zeroth sweep name to first sweep name
$tar_time[0] = $tar_time[1];
$swptime[0] = $swptime[1];
# $s_size[0]  = $s_size[1];  # leave $s_size[0] = 0;
$fx_ang[0]  = $fx_ang[1];
$vol_num[0] = $vol_num[1];
$stype[0]   = $stype[1];
$name[0]    = $name[1];

# dirty hack to replace $name with specific platform name

if ( $my_plat == TRUE ) {
    for($i=0; $i<=$#name; $i++ ) {
	$name[$i] = $this_plat;
    }
}

if ( $show_continuous == TRUE ) {    # special case to list continuous periods of data, then exit.

    print `date`;
    print "\n\tAllowed time_gap = $tgap\n";
    printf("\n\tStart Sweep:\t%s\n\tEnd Sweep:  \t%s\n\n",$sweeps[0],$sweeps[$#sweeps]);

    $first = $sweeps[0];
    for($i=1; $i<=$#sweeps; $i++) {
	$tdiff = $swptime[$i] - $swptime[$i-1];
	if ( abs($tdiff) > $tgap ) {
	    printf("%s\t-> %s\t---> %3d\n", $first, $sweeps[$i-1], $tdiff);
	    $first = $sweeps[$i];
	}
    }
    printf("%s\t-> %s\n\n", $first, $sweeps[$#sweeps]);

    exit;
}

if ( $use_minutes == TRUE ) {	# special case to bundle by fixed clock minute interval

# create yet-another-array; this one assigns each time to a time interval

    for($i=0; $i<=$#sweeps; $i++) {
	push @interv, ($swptime[$i] - ($swptime[$i] % $sec_int));
    }
    
    ($sec,$min,$hr,$day,$mon,$yr)=gmtime($interv[0]);
    $temp = sprintf("%04d%02d%02d_%02d%02d%02d",$yr+1900,$mon+1,$day,$hr,$min,$sec);
    printf("tar -cvf %s%s_%s%s.tar %s",$odir,$name[0],$temp,$fmt,$incl);

    for($i=1; $i<=$#sweeps; $i++) {

	if ( $interv[$i] != $interv[$i-1] ) {
	    ($sec,$min,$hr,$day,$mon,$yr)=gmtime($interv[$i]);
	    $temp = sprintf("%04d%02d%02d_%02d%02d%02d",$yr+1900,$mon+1,$day,$hr,$min,$sec);
	    printf("\n\ntar -cvf %s%s_%s%s.tar %s",$odir,$name[$i],$temp,$fmt,$incl);
	}
	printf(" \\\n\t%s",@sweeps[$i]);

    }
    printf("\n");
    exit;
}


# create yet-another-array; this one is a cumulative size, floored to size intervals

    $tot_size = 0;
    for($i=0; $i<=$#sweeps; $i++) {
	$tot_size += $s_size[$i];
	push @size_int, int($tot_size/$fsize);  # this allows fsize to be more of an average size
    }
    
if ( $use_size_only == TRUE ) {   # special case if size is the sole bundling criteria, then exit

    printf("tar -cvf %s%s_%s%s.tar %s",$odir,$name[0],$tar_time[0],$fmt,$incl);

    for($i=1; $i<=$#sweeps; $i++) {

	if ( $size_int[$i] != $size_int[$i-1] ) {
	    printf("\n\ntar -cvf %s%s_%s%s.tar %s",$odir,$name[$i],$tar_time[$i],$fmt,$incl);
	}
	printf(" \\\n\t%s",@sweeps[$i]);

    }
    printf("\n");
    exit;
}

$arch_start = $swptime[0];

if( $ignore_type == FALSE ) {

    printf("tar -cvf %s%s_%s.%s%s.tar %s",$odir,$name[0],@tar_time[0],@stype[0],$fmt,$incl);

    for($i=1; $i<=$#sweeps; $i++) {

        $vol_chg = ( (@stype[$i] ne @stype[$i-1]) ||
                    (( ($use_vol_num == TRUE ) && ($vol_num[$i] != $vol_num[$i-1]) ) ||
		     ( ($use_vol_num == FALSE) && (@fx_ang[$i] < @fx_ang[$i-1]   ) ) )) ? TRUE : FALSE;

        $tdiff = @swptime[$i] - @swptime[$i-1];

        if ( $tdiff >= $tgap  || (@stype[$i] ne @stype[$i-1]) ||
	     (($vol_chg == TRUE) && ($use_multiv == FALSE)) ||
             (($vol_chg == TRUE) && ((@swptime[$i] - $arch_start) > $archtime ))) {
	    $arch_start = @swptime[$i];
	    printf("\n\ntar -cvf %s%s_%s.%s%s.tar %s",$odir,$name[$i],@tar_time[$i],@stype[$i],$fmt,$incl);
        }
	    printf(" \\\n\t%s",@sweeps[$i]);
    }
    printf("\n");
    exit;
}

if( $ignore_type == TRUE ) {

    printf("tar -cvf %s%s_%s%s.tar %s",$odir,$name[0],@tar_time[0],$fmt,$incl);

    $c_size = 0;
    for($i=1; $i<=$#sweeps; $i++) {

	$c_size += $s_size[$i];
        $vol_chg = ((( ($use_vol_num == TRUE ) && ($vol_num[$i] != $vol_num[$i-1]) ) ||
		     ( ($use_vol_num == FALSE) && (@fx_ang[$i] < @fx_ang[$i-1])) ))  ? TRUE : FALSE;

        $tdiff = @swptime[$i] - @swptime[$i-1];
#	printf("%10d\n",$tdiff);
        if ( $tdiff >= $tgap  || (($vol_chg == TRUE) && ($use_multiv == FALSE)) ||
	     ($use_size == TRUE && ( $c_size >= $fsize )) ||
             (($vol_chg == TRUE) && ((@swptime[$i] - $arch_start) > $archtime ))) {
	    $arch_start = @swptime[$i];
	    $c_size = 0;
	    printf("\n\ntar -cvf %s%s_%s%s.tar %s",$odir,$name[$i],@tar_time[$i],$fmt,$incl);   # don't include stype
        }
	printf(" \\\n\t%s",@sweeps[$i]);
    }
    printf("\n");
    exit;

}
