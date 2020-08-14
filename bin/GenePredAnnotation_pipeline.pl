#!/usr/bin/perl -w
use Getopt::Long qw(GetOptions);

## ------------------------------------------------------- USAGE
my $usage = "\n\tUsage: perl $0 -i <FILE_RUTE_FASTA> -o <OUTPUT FOLDER> <OPTIONS>

		Options:

		**Main analyzes**
		
		-A			Gene annotation (EggNOG mapper and InterproSCAN)
		-P			Gene prediction (Prodigal)
		
		**OTHER**
		-i		Input list of genomes/proteomes
		-o		Output folder		
		-t|--threads	Number of Threads		
		-h		Print help
		-v		version

		**NOTES** 
		1) If you only use annotation step, provide gene prediction in amino acid sequence
		2) To run all analizes please select options -A -P
		3) Always select threads when run annotation steps
		4) To run this pipeline need to install Prodigal, eggNOG mapper and InterproSCAN in your computer\n\n" ;

		
	

## ------------------------------------------------------- DIE if not arguments
if (@ARGV < 1){
	die "\nNot option found!$usage";
}




## ------------------------------------------------------- OPTIONS


my $help;
my $input;
my $threads;
my $annotation;
my $prodigal;
my $outfolder;
my $version;

GetOptions(

	'outfolder|o=s' => \$outfolder,	
	'A' => \$annotation,
	'h' => \$help,
	'input|i=s' => \$input,
	'threads|t=s' => \$threads,
	'P' => \$prodigal,
	'v' => \$version,
	

) or die "$usage";




## ------------------------------------------------------- CHECK OPTIONS

if ($help){
	die "\n\t----------------- You selected help (-h) -----------------\n$usage";
}

if ($version){
	die "GAP v1.0\n";
}

if (!$input) {
	die "\nNot input file. Use -h for help\n\n";
}




## ------------------------------------------------------- MAKING FOLDERS

system "mkdir -p $outfolder/prodigal_outputs";
system "mkdir -p $outfolder/prodigal_outputs/gff_files";
system "mkdir -p $outfolder/prodigal_outputs/aa_fasta";
system "mkdir -p $outfolder/prodigal_outputs/nt_fasta";
system "mkdir -p $outfolder/eggNOG_outputs/raw";
system "mkdir -p $outfolder/eggNOG_outputs/parsed";
system "mkdir -p $outfolder/interproSCAN_outputs/raw";
system "mkdir -p $outfolder/interproSCAN_outputs/parsed";
system "mkdir -p $outfolder/Final_annotation";


## ------------------------------------------------------- SCRIPTS PERL REQUIRED

my $pscript = "GenePredAnnot/bin/parse_interproannotation.py";
my $rscript = "GenePredAnnot/bin/merge_annotations.R";


##                                                             CODE



## --------------------------------------------- 1.- Gene prediction and Annotation----------------------------------------------------- ##

if (defined $prodigal){

	## ------------------------ Prodigal prediction

	my @files = `cat $input`;
	chomp @files;
	
	foreach my $rutes (@files){
		@div_rute = split ("/", $rutes);
        	$name_file1 = $div_rute[-1];
        	@splited_name= split(/\_/,	$name_file1);
        	$outfile_name= join("_", $splited_name[0], $splited_name[1],$splited_name[2]);
		system "prodigal -i $rutes -o $outfolder/prodigal_outputs/gff_files/$outfile_name\.gff -f gff -a $outfolder/prodigal_outputs/aa_fasta/$outfile_name\.faa -d $outfolder/prodigal_outputs/nt_fasta/$outfile_name\.fna\n";
		system "sed -i 's/ #.*//g' $outfolder/prodigal_outputs/aa_fasta/*.faa";
		system "sed -i 's/*//g' $outfolder/prodigal_outputs/aa_fasta/*.faa";
		system "sed -i 's/ #.*//g' $outfolder/prodigal_outputs/nt_fasta/*.fna";
		
	}
	system "find $outfolder/prodigal_outputs/aa_fasta/ -type f | sort > $outfolder/prodigal_outputs/fasta_aa.list";


	if (defined $annotation){

	## ------------------------ interpro scan

		my @files2 = `cat $outfolder/prodigal_outputs/fasta_aa.list`;
		chomp @files2;

		foreach my $rutes2 (@files2){
			@div_rute2 = split ("/", $rutes2);
        		$name_file2 = $div_rute2[-1];
        		@splited_name2= split(/\_/,	$name_file2);
			@splited_final= split(/\./, $splited_name2[2]);
        		$outfile_name2= join("_", $splited_name2[0], $splited_name2[1],$splited_final[0]);
			system "interproscan.sh -i $rutes2  -goterms -pa -f tsv -o $outfolder/interproSCAN_outputs/raw/$outfile_name2\.tsv -cpu $threads\n";
			
		}

	## ---------------------- eggNOG annotation
		my @files3 =`cat $outfolder/prodigal_outputs/fasta_aa.list`;
		chomp @files3;
	
		foreach my $rutes3 (@files3){
			@div_rute3 = split ("/", $rutes3);
        		$name_file3 = $div_rute3[-1];
        		@splited_name3= split(/\_/,	$name_file3);
			@splited_final2= split(/\./, $splited_name3[2]);
        		$outfile_name3= join("_", $splited_name3[0], $splited_name3[1],$splited_final2[0]);
			system "emapper.py -i $rutes3  --output $outfolder/eggNOG_outputs/raw/$outfile_name3 -m diamond --cpu $threads\n";

			
		}
	##------------------- parse interproscan outputs
		system "sed -i 1i\'ID\tSequence_MD5\tLength\tAnalysis\tAnAccess\tDescription\tStart\tStop\tScore\tStatus\tDate\tInterPro_ann\tInterpro_description\tGO\tPathways\' $outfolder/interproSCAN_outputs/raw/*";	
		system "find $outfolder/interproSCAN_outputs/raw/ -type f | sort > $outfolder/interproSCAN_outputs/interpro.list";
		my @interprOut = `cat $outfolder/interproSCAN_outputs/interpro.list`;
		chomp @interprOut;
		
		foreach my $interIn (@interprOut){
			@div= split ("/", $interIn);
			$fn= $div[-1];
        		@sname= split(/\_/,	$fn);
			@sfinal= split(/\./, $sname[2]);
			$outname= join("_", $sname[0],$sname[1],$sfinal[0]);
			system "python3 $pscript $interIn $outfolder/interproSCAN_outputs/parsed/$outname\_parsed.tsv\n";
	
		}
		
	##----------------- parse eggNOG outputs
		system "mkdir -p $outfolder/eggNOG_outputs/raw/annotations $outfolder/eggNOG_outputs/raw/seed_orthologs"; 
		system "mv $outfolder/eggNOG_outputs/raw/*.annotations $outfolder/eggNOG_outputs/raw/annotations";
		system "mv $outfolder/eggNOG_outputs/raw/*.seed_orthologs $outfolder/eggNOG_outputs/raw/seed_orthologs";
		system "find $outfolder/eggNOG_outputs/raw/annotations -type f | sort > $outfolder/eggNOG_outputs/eggNOG.list";
		my @eggOut =`cat $outfolder/eggNOG_outputs/eggNOG.list`;
		chomp @eggOut;
	
		foreach my $eggIn (@eggOut){
			@div2= split ("/", $eggIn);
			$fn2= $div2[-1];
			@sname2= split (/\_/, $fn2);
			@sfinal2= split (/\./, $sname2[2]);
			$outname2= join("_", $sname2[0],$sname2[1],$sfinal2[0]);
   			system "sed '/^#/ d' $eggIn > $outfolder/eggNOG_outputs/parsed/$outname2\.tmp\n";
		}
		
		system "sed -i 1i\'gene_id\tseed_eggNOG_ortholog\tseed_ortholog_evalue\tseed_ortholog_score\tbest_tax_level\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tTax.scope\teggNOG.Ogs\tBest.Ogs\tCOGs\teggNOG.HMM.Description\' $outfolder/eggNOG_outputs/parsed/*\n" ;

		

	

	##--------------------- join eggnog interpro annotation		
		system "find $outfolder/interproSCAN_outputs/parsed/ -type f | sort > $outfolder/interproSCAN_outputs/interpro_parsed.list";
		system "find $outfolder/eggNOG_outputs/parsed/	 -type f | sort > $outfolder/eggNOG_outputs/eggnog_parsed.list";				
		system "paste $outfolder/eggNOG_outputs/eggnog_parsed.list $outfolder/interproSCAN_outputs/interpro_parsed.list | sed 's/ /\t/g' > $outfolder/Final_annotation/final.list";		
		my @inlist = `cat $outfolder/Final_annotation/final.list`;
		chomp @inlist;

		foreach my $finalist (@inlist){
		($rute_file1, $rute_file2) = split ("\t", $finalist);
		system "Rscript --vanilla $rscript $rute_file1 $rute_file2 $outfolder/Final_annotation\n";
		}

	}
}
	







#------------------------------------------------------2. Gene function annotation ---------------------------------------------------------
if(!$prodigal){	



	##----------interproscan annotation

my @files4 = `cat $input`;
	chomp @files4;

	foreach my $rutes4 (@files4){
		@div_rute4 = split ("/", $rutes4);
        	$name_file4 = $div_rute4[-1];
        	@splited_name4= split(/\_/,	$name_file4);
		@splited_final4= split(/\./, $splited_name4[2]);
        	$outfile_name4= join("_", $splited_name4[0], $splited_name4[1],$splited_final4[0]);
		system "interproscan.sh -i $rutes4  -goterms -pa -f tsv -o $outfolder/interproSCAN_outputs/raw/$outfile_name4\.tsv -cpu $threads\n";

	}

	## ---------------------- eggNOG annotation
	my @files5 =`cat $input`;
	chomp @files5;
	
	foreach my $rutes5 (@files5){
		@div_rute5 = split ("/", $rutes5);
        	$name_file5 = $div_rute5[-1];
        	@splited_name5= split(/\_/,	$name_file5);
		@splited_final5= split(/\./, $splited_name5[2]);
        	$outfile_name5= join("_", $splited_name5[0], $splited_name5[1],$splited_final5[0]);
		system "emapper.py -i $rutes5  --output $outfolder/eggNOG_outputs/raw/$outfile_name5\.tsv -m diamond --cpu $threads\n";

	}
	##---------------------- parse interproScan output
	system "sed -i 1i\'ID\tSequence_MD5\tLength\tAnalysis\tAnAccess\tDescription\tStart\tStop\tScore\tStatus\tDate\tInterPro_ann\tInterpro_description\tGO\tPathways\' $outfolder/interproSCAN_outputs/raw/*";
	system "find $outfolder/interproSCAN_outputs/raw/ -type f | sort > $outfolder/interproSCAN_outputs/interpro.list";
	my @interprOut = `cat $outfolder/interproSCAN_outputs/interpro.list`;
	chomp @interprOut;
		
	foreach my $interIn (@interprOut){
		@div= split ("/", $interIn);
		$fn= $div[-1];
        	@sname= split(/\_/,	$fn);
		@sfinal= split(/\./, $sname[2]);
		$outname= join("_", $sname[0],$sname[1],$sfinal[0]);	
		system "python3 $pscript $interIn $outfolder/interproSCAN_outputs/parsed/$outname\_parsed.tsv\n";
		
		
	}
		

	##-------------------- parse eggnog output
		system "mkdir -p $outfolder/eggNOG_outputs/raw/annotations $outfolder/eggNOG_outputs/raw/seed_orthologs"; 
		system "mv $outfolder/eggNOG_outputs/raw/*.annotations $outfolder/eggNOG_outputs/raw/annotations";
		system "mv $outfolder/eggNOG_outputs/raw/*.seed_orthologs $outfolder/eggNOG_outputs/raw/seed_orthologs";
		system "find $outfolder/eggNOG_outputs/raw/annotations -type f | sort > $outfolder/eggNOG_outputs/eggNOG.list";
		my @eggOut =`cat $outfolder/eggNOG_outputs/eggNOG.list`;
		chomp @eggOut;
	
		foreach my $eggIn (@eggOut){
			@div2= split ("/", $eggIn);
			$fn2= $div2[-1];
			@sname2= split (/\_/, $fn2);
			@sfinal2= split (/\./, $sname2[2]);
			$outname2= join("_", $sname2[0],$sname2[0],$sfinal2[0]);
   			system "sed '/^#/ d' $eggIn > $outfolder/eggNOG_outputs/parsed/$outname2\.tmp \n";

		}
		
		system "sed -i 1i\'gene_id\tseed_eggNOG_ortholog\tseed_ortholog_evalue\tseed_ortholog_score\tbest_tax_level\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tTax.scope\teggNOG.Ogs\tBest.Ogs\tCOGs\teggNOG.Description\' $outfolder/eggNOG_outputs/parsed/*\n" ;

		
		


	##--------------------- join eggnog interpro annotation
	system "find $outfolder/interproSCAN_outputs/parsed/ -type f | sort > $outfolder/interproSCAN_outputs/interpro_parsed.list";
	system "find $outfolder/eggNOG_outputs/parsed/	 -type f | sort > $outfolder/eggNOG_outputs/eggnog_parsed.list";	
	system "paste $outfolder/eggNOG_outputs/eggnogparsed.list $outfolder/interproSCAN_outputs/interpro_parsed.list | sed 's/ /\t/g' > $outfolder/Final_annotation/final.list";			
	my @inlist = `cat $outfolder/Final_annotation/final.list`;
	chomp @inlist;

	foreach my $finalist (@inlist){
	($rute_file1, $rute_file2) = split ("\t", $finalist);
	system "Rscript --vanilla $rscript $rute_file1 $rute_file2 $outfolder/Final_annotation\n";
	}

}
