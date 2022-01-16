#!/usr/bin/perl
##
## creates a child process per chain
##usage perl forkFilter.pl #forks
use warnings;
use strict;
use Parallel::ForkManager;
#
my $max = shift(@ARGV); #get number of cores to use at one time
my $pm = Parallel::ForkManager->new($max);

my $dir = '/uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/SV_files_new/bamfiles/entropyrunnew/';
my $odir = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/nanoporereads/SV_files_new/bamfiles/entropyrunnew/outfiles/";

CHAINS:

foreach my $j (0..4){
                 system"sleep 3\n";
                 $pm->start and next CHAINS; ##fork;
				print "/uufs/chpc.utah.edu/common/home/u6000989/bin/bgc_v1.05b -A $dir"."LI.combSNPs.gl -B $dir"."LM.combSNPs.gl -Z $dir"."DBS.combSNPs.gl -a $dir"."LI.LKfilter4new.gl -b $dir"."LM.LKfilter4new.gl -h $dir"."DBS.LKfilter4new.gl -F $odir"."bgcout_pct3_filter4_$j -O 0 -x 200000 -n 20000 -t 5 -p 1 -q 1 -N 1 -m 0 -d 1 -s 0 -I 1 -T 1.8 -R 1";
                          	system "/uufs/chpc.utah.edu/common/home/u6000989/bin/bgc_v1.05b -A $dir"."LI.combSNPs.gl -B $dir"."LM.combSNPs.gl -Z $dir"."DBS.combSNPs.gl -a $dir"."LI.LKfilter4new.gl -b $dir"."LM.LKfilter4new.gl -h $dir"."DBS.LKfilter4new.gl -F $odir"."bgcout_pct3_filter4_$j -O 0 -x 200000 -n 20000 -t 5 -p 1 -q 1 -N 1 -m 0 -d 1 -s 0 -I 1 -T 1.8 -R 1";
                $pm->finish;
}
$pm->wait_all_children;

