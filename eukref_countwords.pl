#!/usr/bin/perl

# script to count words in GenBank records from sequences retrieved by eukref_gbretrival.py
# Daniel Richter, 21 July 2015

use warnings;
use strict;

my $DIRECTORY_DELIMITER = "/";

my $RESULTS_FASTA_FILE = "temp_results_silva_uchime.fas";

opendir(DIR, ".") || die "could not open current directory for read";

my %words_by_cycle = ();
my %word_counts = ();

my $max_cycle = 0;

foreach my $subdirectory (readdir(DIR))
{
    if ($subdirectory =~ /^cycle_(\d+)/)
    {
	my $cycle = $1;

	if ($cycle > $max_cycle)
	{
	    $max_cycle = $cycle;
	}

	my $fasta_file_to_read = $subdirectory . $DIRECTORY_DELIMITER . $RESULTS_FASTA_FILE;

	open(INPUT, $fasta_file_to_read) || die "could not open '$fasta_file_to_read' for read";
	
	while (<INPUT>)
	{
	    if (/^>/)
	    {
		chomp;

		my (undef, undef, undef, undef, $description) = split(/\|/);

		if ($description =~ /^\s+(.*)/)
		{
		    $description = $1;
		}

		$description =~ s/,//g;
		$description = lc $description;

#		print STDERR $description . "\n";

		foreach my $word (split(/\s+/, $description))
		{
		    if (not defined $word_counts{$word})
		    {
			$word_counts{$word} = 1;
		    }
		    else
		    {
			$word_counts{$word}++;
		    }

		    if (not defined $words_by_cycle{$word}->{$cycle})
		    {
			$words_by_cycle{$word}->{$cycle} = 1;
		    }
		    else
		    {
			$words_by_cycle{$word}->{$cycle}++;
		    }
		}
	    }
	}

	close INPUT || die "could not close '$fasta_file_to_read' after read";
    }
}

closedir DIR;

print "word";

for (my $cycle = 1; $cycle <= $max_cycle; $cycle++)
{
    print "\t" . $cycle;
}

print "\n";

foreach my $word (sort { $word_counts{$b} <=> $word_counts{$a} } keys %word_counts)
{
    print $word;

    for (my $cycle = 1; $cycle <= $max_cycle; $cycle++)
    {
	if (defined $words_by_cycle{$word}->{$cycle})
	{
	    print "\t" . $words_by_cycle{$word}->{$cycle};
	}
	else
	{
	    print "\t0";
	}
    }

    print "\n";
}
