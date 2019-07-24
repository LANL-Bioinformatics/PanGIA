#!/usr/bin/env perl
use Cwd;
use strict;

my $bin_dir = Cwd::realpath($0);
$bin_dir =~ s/pangia-vis\.pl//;
my $pvis_dir = "$bin_dir/pangia-vis";
my $tsvfile = Cwd::realpath($ARGV[0]);

# get project prefix
if( @ARGV > 0 )
{
    my ($prefix) = $tsvfile =~ /^(.*)\.report\.tsv/;

    if( -e $tsvfile && $tsvfile=~/\.report\.tsv/ ){
        if( -d "${prefix}_tmp" ){
            if( -d $prefix ){
                print STDERR "[WARNING] Genome coverage directory existed: $prefix/.\n";
            }
            else{
                print STDERR "[INFO] Generating genome coverage data...\n";
                my $ecode = system("$pvis_dir/scripts/depth_scale_down.sh ${prefix}_tmp/merged_sam $prefix");
                if( $ecode ){
                    print STDERR "[ERROR] Failed to calculate genome coverage data.\n";
                }
                else{
                    print STDERR "[INFO] Done.\n";
                }
            }
        }
        else{
            print STDERR "[WARNING] ${prefix}_tmp directory not found. Genome coverage map is unavailable.\n";
        }
    }
    else{
        die( "[ERROR] Input PanGIA result not found: $tsvfile\n" );
    }
    print STDERR "[INFO] Opening PanGIA-VIS application on http://localhost:5006/pangia-vis\n";
    `bokeh serve $pvis_dir --show --args $tsvfile 1>&2`;
}
else{
    die("USAGE: ./pangia-vis.pl path_to/pangia.report.tsv\n");
}
