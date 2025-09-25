#!/usr/bin/env perl

BEGIN
{
    push(@INC, (getpwnam('drichter'))[7] . "/lib/perl");
}

use warnings;
use strict;

use GetOptions;

my $HELP_OPTION = 'help';
my $QUERY_OPTION = 'query';
my $DATABASE_OPTION = 'database';
my $OUTPUT_OPTION = 'output';

my %OPTION_TYPES = ($HELP_OPTION => '',
		    $QUERY_OPTION => '=s',
		    $DATABASE_OPTION => '=s',
		    $OUTPUT_OPTION => '=s');

sub main 
{
    my %args = &GetOptions::get_options(\%OPTION_TYPES);

    if ($args{$HELP_OPTION} || (not ($args{$QUERY_OPTION} && $args{$DATABASE_OPTION} && $args{$OUTPUT_OPTION})))
    {
	&help();
	exit;
    }

    foreach my $directory ($args{$QUERY_OPTION}, $args{$DATABASE_OPTION}, $args{$OUTPUT_OPTION})
    {
	if (not -d $directory)
	{
	    die "'$directory' is not a directory";
	}
    }

    opendir(QUERY, $args{$QUERY_OPTION}) || die "could not open directory '$args{$QUERY_OPTION}' for read";
    opendir(DATABASE, $args{$DATABASE_OPTION}) || die "could not open directory '$args{$DATABASE_OPTION}' for read";

    my @query_files = readdir(QUERY);
    my @database_files = readdir(DATABASE);

    closedir QUERY;
    closedir DATABASE;
    
    foreach my $query_file (sort @query_files)
    {
	if ($query_file !~ /\.fasta$/)
	{
	    next;
	}

	foreach my $database_file (sort @database_files)
	{
	    if ($database_file !~ /\.fna$/)
	    {
		next;
	    }
	    
	    print STDOUT "--- $query_file vs. $database_file ---\n";
	    
	    my $output_path = $args{$OUTPUT_OPTION} . "/" .
		$query_file . "_" . $database_file . "_blastout_decontam_refG_out6.txt";

	    if (-e $output_path)
	    {
		print STDOUT "skipping, output file $output_path exists\n";

		next;
	    }

	    my $query_path = $args{$QUERY_OPTION} . "/" . $query_file;
	    my $database_path = $args{$DATABASE_OPTION} . "/" . $database_file;
	    
	    my $blastn_command = "blastn -query $query_path -db $database_path -max_target_seqs 1 -outfmt 6 -num_threads 8 " .
		"-out $output_path";

	    print STDOUT $blastn_command . "\n";

	    `$blastn_command`;
	}
    }
}

sub help
{
    my $HELP = <<HELP;
Syntax: $0 -$QUERY_OPTION <dir> -$DATABASE_OPTION <dir> -$OUTPUT_OPTION <dir>

Run blastn for all of the files in the query directory against all of
the files in the database directory.

    -$HELP_OPTION : print this message
    -$QUERY_OPTION : directory containing query .fasta files
    -$DATABASE_OPTION : directory containing query .fna files
    -$OUTPUT_OPTION : output directory

HELP

    print STDERR $HELP;

}

&main();
