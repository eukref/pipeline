#!/opt/local/bin/perl -w

#Quick and dirty perl script to get certain metainformation from genebank accession numbers.

# Mallo, D. 2012; del Campo, J. 2015

use strict;
use warnings;
use Data::Dumper;
use LWP::Simple ;
use Time::HiRes qw(usleep nanosleep gettimeofday);
#use Bio::DB::GenBank;

#Objects
#my $gb=Bio::DB::GenBank->new(-complexity=>1);
my $seq_ob;
my $file;
my $filename;
my @ans;
my ($feat,$val);
my $outfile;

#Numerical values
my $max_it=10;
my $i=0;
my $downloaded=0;

#Strings
my ($cln,$org,$mol,$cou,$iso,$env,$hst);


#$gb->retrieval_type('io_string'); #############VERY IMPORTANT!!!!: The default value 'pipeline' uses fork and pipes and these generate huge problems being inside the eval!!!!!

if (scalar(@ARGV)!=1 || !(-f $ARGV[0]))
{
	die "Error, you should run this script as: $0 accession_number_list_file\n";
}

$filename=$ARGV[0];

open ($file,$filename);
@ans=<$file>;
close $file;

open ($outfile,">$filename\_out\.csv");
my $tax ; 
foreach my $an (@ans)
{
chomp $an ;
my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$an&retmode=xml" ;
 
    my $tax ; 
    my $content = get($url);
#    print $content,"\n" ;
    my $featrefs = GetFeaturesQualifiers($content );
    
    $org = GetTag($content , "GBSeq_organism") ; 
    if( defined($featrefs->{"clone"})){
	$cln = $featrefs->{"clone"} ;
    } else {
	$cln = "n/a" ; 
    }
    $mol = GetTag($content , "GBSeq_moltype") ; 
    $cou = GetTag($content , "GBSeq_country") ; 
    $tax = GetTag($content , "GBSeq_taxonomy") ; 

    if( defined($featrefs->{"isolation_source"})){
	$iso = $featrefs->{"isolation_source"} ;
    } else {
	$iso = "n/a" ; 
    }
    if( defined($featrefs->{"environmental_sample"})){
	$env = "TRUE" ;
    } else {
	$env = "n/a" ; 
    }
    if( defined($featrefs->{"host"})){
	$hst = $featrefs->{"host"} ;
    } else {
	$hst = "n/a" ; 
    }


    print $outfile "$an\t$cln\t$org\t$tax\t$env\t$mol\t$cou\t$hst\t$iso\n";
    print "$an\t$cln\t$org\t$tax\t$env\t$mol\t$cou\t$hst\t$iso\n";


usleep(300) ;

}

sub GetFeaturesQualifiers{
    my( $content) = @_ ;
    my %hashref = () ;
    while( $content =~ m/<GBQualifier>(.+?)<\/GBQualifier>/msg){
	my $tname = GetTag($1,"GBQualifier_name") ; 
	my $tvalue = GetTag($1,"GBQualifier_value") ; 
     #  print "$tname -> $tvalue\n" ;
	$hashref{$tname} = $tvalue ;  
    }
    return \%hashref ;
}
sub GetInTag{
    my( $content , $tag , $what ) = @_ ;
    
}
sub GetTag{
    my( $content,$tag) = @_ ;
    if($content =~ m/<$tag>(.+)<\/$tag>/){
	return $1 ;
    } else {
	return "n/a" ;
    }
}

exit(1) ; 
foreach my $an (@ans)
{

	$i=0;
	$downloaded=0;
	chomp $an;
	print "\nAccesion number to get: $an\n";
	
	#Genbank search (NCBI)
	while($downloaded==0)
	{
		eval
		{
#			$seq_ob = $gb->get_Seq_by_acc($an);
		};
		if (!($@)) 
		{
			print "Downloaded correctly\n";
      		$downloaded=1;
   		}
		else
		{
			$seq_ob=undef;
			if ($i==$max_it)
			{
				die "Max number of download attempts reached, net error\n";
			}
			print "Retrying the download, attempt $i of $max_it\n";
			++$i;
			sleep 2;
		}
		
	}
	#print "Toda la secuencia\n";
	#print Dumper $seq_ob;
	($cln,$org,$mol,$cou,$iso,$env,$hst)=(undef,undef,undef,undef,undef,undef,undef,undef);

	FEAT: for $feat($seq_ob->get_SeqFeatures) 
	{
		#print "Nueva feature\n";
		#print Dumper $feat;
		if($feat->has_tag("clone"))  
    	{
        	for $val($feat->get_tag_values("clone")) 
        	{
        		$cln.="$val";
        	}
        }
    	if($feat->has_tag("organism"))  
    	{
        	for $val($feat->get_tag_values("organism")) 
        	{
        		$org.="$val";
        	}
        }

        if($feat->has_tag("mol_type"))  
    	{
        	for my $val($feat->get_tag_values("mol_type")) 
        	{
        		$mol.="$val";
        	}
        }

        if($feat->has_tag("country"))  
    	{
        	for my $val($feat->get_tag_values("country")) 
        	{
        		$cou.="$val";
        	}
        }

        if($feat->has_tag("isolation_source"))  
    	{
        	for my $val($feat->get_tag_values("isolation_source")) 
        	{
        		$iso.="$val";
        	}
        }

        if($feat->has_tag("environmental_sample"))  
    	{
        	for my $val($feat->get_tag_values("environmental_sample")) 
        	{
        		$env.="1";
        	}
        }
        if($feat->has_tag("host"))  
    	{
        	for $val($feat->get_tag_values("host")) 
        	{
        		$hst.="$val";
        	}
        }
  
    }
    
    !$an and $an='n/a';
    !$cln and $cln='n/a';
    !$org and $org='n/a';
    !$mol and $mol='n/a';
    !$cou and $cou='n/a';
    !$iso and $iso='n/a';
    !$env and $env='0';
    !$hst and $hst='n/a';

    print $outfile "$an;$cln;$org;$tax;$env;$mol;$cou;$hst;$iso\n";
    print "$an;$cln;$org;$tax;$env;$mol;$cou;$hst;$iso\n";
	#exit; #DEBUG
}
close $outfile;
exit;
# 
# ACCESSION
# SOURCE
# ORGANISM 
# FEATURES
# /organism
# /mol_type
# /country
# /isolation_source
# /environmental_sample
# /host
