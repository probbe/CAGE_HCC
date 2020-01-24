#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
use Time::HiRes qw( time );
use Getopt::Long;
use List::Util qw (sum shuffle min max);
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#
#	Description
#		This is a perl script convert bam files into CTSS bed.
#
#	Input
#		--bamPath=			file path[compulsory]; path of the BED file contains the location of the phase1 raw DPI clusters; 
#		--minLen=			int [18]; minimum length to be taken;
#		--maxLen=			int [100]; max length to be taken;
#		--min_MAPQ=			int [20]; minimum min_MAPQ to take;
#		--longShortFormat=	long or short [short]; long = expand the read 
#		--exclude_flag=		string [100]; comma delimited sam bits
#		--outDir=			directory path ['./output/']; output directory;
#
#	v0.1
#	[January 17, 2015 20:24 PM] debut;
#
#	v0.2
#	[June 16, 2015 12:55 PM] added minLen and maxLen option
#
#	v0.3
#	[April 28, 2017 11:35 AM] added outputPrefix
#
#	v0.4
#	[October 4, 2017 18:37 PM] added exclude_flag option
#
#	v0.5
#	[November 3, 2017 16:16 PM] added longShortFormat option
#
#	v0.6
#	[June 4, 2018 17:07 PM] added min_MAPQ option for Marina's project
#	
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	/home/hon-chun/resources/perlScript/FANTOM/bamToCTSS/v0.6/bamToCTSS_v0.6.pl
#	--bamPath=/home/m.lizio/transposable_elements_fantom5/bam/293SLAM%20rinderpest%20infection%2c%2000hr%2c%20biol_rep1.CNhs14406.13541-145H4.hg19.nobarcode.bam
#	--minLen=0
#	--maxLen=100
#	--min_MAPQ=20
#	--longShortFormat=short
#	--exclude_flag=512
#	--outputPrefix=MAPQ_20
#	--outDir=/home/hon-chun/resources/perlScript/FANTOM/bamToCTSS/demo_output_2
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $scriptDirPath = dirname(rel2abs($0));
my $scriptAbsPath = abs_path($0);
my ($curntTimeStamp) = &timeStamp();#->628
my $globalTmpLogPath = "$scriptDirPath/00_screen_status.$curntTimeStamp.log.txt";
my $ARGVStr = join "\n", (&currentTime(), $scriptAbsPath, @ARGV);#->276
my $globalReadmeHsh_ref = {};
open TMPLOG, ">", $globalTmpLogPath;
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#
#<section ID="startingTasks" num="0">
&printCMDLogOrFinishMessage("CMDLog");#->421
my ($bamPath, $minLen, $maxLen, $exclude_flag, $longShortFormat, $min_MAPQ, $outputPrefix, $outDir) = &readParameters();#->569
&checkSamtoolsVersion();#->178
&checkBedtoolsVersion();#->158
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#
#<section ID="defineHardCodedParam" num="1">
my $paramTag = "$outputPrefix";
my $bedtools_bin = "/home/hon-chun/resources/bin/bedtools/bedtools2-2.20.1/bin//bedtools";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag"; push @mkDirAry, $resultDir;
my $resultBedDir = "$resultDir/bed/"; push @mkDirAry, $resultBedDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
my $resultScriptDir = "$resultDir/script/"; push @mkDirAry, $resultScriptDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_convert
#
#<section ID="convert" num="3">
my ($chromNameHsh_ref) = &generateChromNameConversionHsh();#->294
&convertBamToCTSS($bamPath, $paramTag, $minLen, $maxLen, $resultBedDir, $chromNameHsh_ref, $exclude_flag, $bedtools_bin, $min_MAPQ, $longShortFormat, $resultLogDir);#->198
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_finishingTasks
#
#<section ID="finishingTasks" num="4">
&logCalledCMDAndScript($ARGVStr, $resultScriptDir, $scriptAbsPath);#->396
&printOutputFileListAndReadme($ARGVStr, $paramTag, $outDir);#->454
&printCMDLogOrFinishMessage("finishMessage");#->421
&copyGlobalStatusLog($globalTmpLogPath, $resultLogDir, $scriptDirPath);#->256
close TMPLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	checkTools [n=2]:
#		checkBedtoolsVersion, checkSamtoolsVersion
#
#	general [n=7]:
#		checkBedtoolsVersion, checkSamtoolsVersion, currentTime
#		logCalledCMDAndScript, printCMDLogOrFinishMessage, readParameters
#		timeStamp
#
#	log [n=2]:
#		copyGlobalStatusLog, reportAndLogStatus
#
#	output [n=1]:
#		printOutputFileListAndReadme
#
#	specific [n=1]:
#		convertBamToCTSS
#
#	time [n=1]:
#		timeStamp
#
#	unassigned [n=1]:
#		generateChromNameConversionHsh
#
#====================================================================================================================================================#

sub checkBedtoolsVersion {
#....................................................................................................................................................#
#	subroutineCategory: general, checkTools
#	dependOnSub: reportAndLogStatus|606
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|73
#	secondaryAppearInSection: >none
#	input: none
#	output: none
#	toCall: &checkBedtoolsVersion();
#	calledInLine: 79
#....................................................................................................................................................#

	my $stdOut = `bedtools --version 2>&1`;
	if ($stdOut =~ m/bedtools v(\S+)/) {
		&reportAndLogStatus("Checking: bedtools version: $1", 0, "\n");#->606
	} else {
		die "bedtools not installed properly. Quitting.\n";
	}
}
sub checkSamtoolsVersion {
#....................................................................................................................................................#
#	subroutineCategory: general, checkTools
#	dependOnSub: reportAndLogStatus|606
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|73
#	secondaryAppearInSection: >none
#	input: none
#	output: none
#	toCall: &checkSamtoolsVersion();
#	calledInLine: 78
#....................................................................................................................................................#

	my $samtoolsStdout = `samtools 2>&1`;
	if ($samtoolsStdout =~ m/\s+(Version: \S+)\s+/) {
		&reportAndLogStatus("Checking: samtools: $1", 0, "\n");#->606
	} else {
		die "samtools not installed properly. Quitting.\n";
	}
}
sub convertBamToCTSS {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportAndLogStatus|606
#	appearInSub: >none
#	primaryAppearInSection: 3_convert|105
#	secondaryAppearInSection: >none
#	input: $bamPath, $chromNameHsh_ref, $maxLen, $minLen, $outputPrefix, $resultBedDir
#	output: 
#	toCall: &convertBamToCTSS($bamPath, $outputPrefix, $minLen, $maxLen, $resultBedDir, $chromNameHsh_ref);
#	calledInLine: 109
#....................................................................................................................................................#
	my ($bamPath, $paramTag, $minLen, $maxLen, $resultBedDir, $chromNameHsh_ref, $exclude_flag, $bedtools_bin, $min_MAPQ, $longShortFormat, $resultLogDir) = @_;
	
	my $CTSSBedPath = "$resultBedDir/$paramTag.$longShortFormat.ctss.bed.gz";
	my $logPath = "$resultLogDir/$paramTag.$longShortFormat.read_count.txt";
	&reportAndLogStatus("converting $paramTag bam to ctss bed", 10, "\n");#->606

	#---[1/17/15 21:33] -view -F 0x100 to include primary alignments only
	#open (BAMIN, "samtools view -bF 0x100 $bamPath | bedtools bamtobed -bed12 -split -i stdin |");
	my $num_proc = 0;
	my $F_param = '';
	if (defined $exclude_flag) {
		my @flagAry = ();
		foreach my $flag (split /,/, $exclude_flag) {
			push @flagAry, "-F $flag";
		}
		$F_param = join " ", @flagAry;
	}
	
	#---[6/18/15 14:23] get all
	open (BAMIN, "samtools view $F_param -b $bamPath | bedtools bamtobed -bed12 -split -i stdin |");
	if ($longShortFormat eq 'short') {
		open (CTSSOUT, "| sort -k1,1 -k2,2n | $bedtools_bin merge -d -1 -s -c 4,5,6 -o distinct,count,distinct -i stdin | gzip -c >$CTSSBedPath");
	} elsif ($longShortFormat eq 'long') {
		open (CTSSOUT, "| sort -k1,1 -k2,2n | gzip -c >$CTSSBedPath");
	} else {
		die "--longShortFormat= has to be long or short\n";
	}
	
	my $total_read_num = 0;
	my $passed_read_num = 0;
	
	while (my $line = <BAMIN>) {
		chomp $line;
		my ($chrom, $bedStart, $bedEnd, undef, $MAPQ, $strand, undef, undef, undef, undef, $blockSizes, undef) = split /\t+/, $line;
		$total_read_num++;
		
		next if $MAPQ < $min_MAPQ;

		if (exists $chromNameHsh_ref->{$chrom}) {
			$chrom = $chromNameHsh_ref->{$chrom};
		}
		
		my ($end5Start, $end5End);
		my $length = sum(split/,/, $blockSizes);

		if ($length >= $minLen and $length <= $maxLen) {
			$passed_read_num++;
			if ($strand eq '+') {
				$end5Start = $bedStart;
				$end5End = $end5Start+1;

			} elsif ($strand eq '-') {
				$end5Start = $bedEnd-1;
				$end5End = $bedEnd;
			}

			#my $CTSSID = $chrom.":".$end5Start."..".$end5End.",".$strand;
			my $CTSSID = 'ctss';
			print CTSSOUT join "", (join "\t", ($chrom, $end5Start, $end5End, $CTSSID, '1', $strand)), "\n";
		}
	}
	
	&reportAndLogStatus("total_read_num = $total_read_num", 10, "\n");
	&reportAndLogStatus("passed_read_num = $passed_read_num", 10, "\n");

	open (LOG, ">$logPath");
	print LOG "passed_read_num = $passed_read_num\n";
	print LOG "total_read_num  = $total_read_num\n";
	close LOG;
	
	close BAMIN;
	close CTSSOUT;

	return ();
}
sub copyGlobalStatusLog {
#....................................................................................................................................................#
#	subroutineCategory: log
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|114
#	secondaryAppearInSection: >none
#	input: $globalTmpLogPath, $resultLogDir, $scriptDirPath
#	output: 
#	toCall: &copyGlobalStatusLog($globalTmpLogPath, $resultLogDir, $scriptDirPath);
#	calledInLine: 120
#....................................................................................................................................................#
	my ($globalTmpLogPath, $resultLogDir, $scriptDirPath) = @_;
	
	system "cp -f $globalTmpLogPath $resultLogDir/00_screen_status.this_run.log.txt";
	system "cp -f $globalTmpLogPath $scriptDirPath/00_screen_status.most_recent.log.txt";
	system "rm -f $globalTmpLogPath";
	
	return ();
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|421, reportAndLogStatus|606
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|73, 4_finishingTasks|114
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 62, 441, 444, 449, 622, 623
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateChromNameConversionHsh {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_convert|105
#	secondaryAppearInSection: >none
#	input: none
#	output: $chromNameHsh_ref
#	toCall: my ($chromNameHsh_ref) = &generateChromNameConversionHsh();
#	calledInLine: 108
#....................................................................................................................................................#
	
	my $chromNameHsh_ref = {
		'GL000192.1' => 'chr1_gl000192_random',
		'GL000225.1' => 'chrUn_gl000225',
		'GL000194.1' => 'chr4_gl000194_random',
		'GL000193.1' => 'chr4_gl000193_random',
		'GL000200.1' => 'chr9_gl000200_random',
		'GL000222.1' => 'chrUn_gl000222',
		'GL000212.1' => 'chrUn_gl000212',
		'GL000195.1' => 'chr7_gl000195_random',
		'GL000223.1' => 'chrUn_gl000223',
		'GL000224.1' => 'chrUn_gl000224',
		'GL000219.1' => 'chrUn_gl000219',
		'GL000205.1' => 'chr17_gl000205_random',
		'GL000215.1' => 'chrUn_gl000215',
		'GL000216.1' => 'chrUn_gl000216',
		'GL000217.1' => 'chrUn_gl000217',
		'GL000199.1' => 'chr9_gl000199_random',
		'GL000211.1' => 'chrUn_gl000211',
		'GL000213.1' => 'chrUn_gl000213',
		'GL000220.1' => 'chrUn_gl000220',
		'GL000218.1' => 'chrUn_gl000218',
		'GL000209.1' => 'chr19_gl000209_random',
		'GL000221.1' => 'chrUn_gl000221',
		'GL000214.1' => 'chrUn_gl000214',
		'GL000228.1' => 'chrUn_gl000228',
		'GL000227.1' => 'chrUn_gl000227',
		'GL000191.1' => 'chr1_gl000191_random',
		'GL000208.1' => 'chr19_gl000208_random',
		'GL000198.1' => 'chr9_gl000198_random',
		'GL000204.1' => 'chr17_gl000204_random',
		'GL000233.1' => 'chrUn_gl000233',
		'GL000237.1' => 'chrUn_gl000237',
		'GL000230.1' => 'chrUn_gl000230',
		'GL000242.1' => 'chrUn_gl000242',
		'GL000243.1' => 'chrUn_gl000243',
		'GL000241.1' => 'chrUn_gl000241',
		'GL000236.1' => 'chrUn_gl000236',
		'GL000240.1' => 'chrUn_gl000240',
		'GL000206.1' => 'chr17_gl000206_random',
		'GL000232.1' => 'chrUn_gl000232',
		'GL000234.1' => 'chrUn_gl000234',
		'GL000202.1' => 'chr11_gl000202_random',
		'GL000238.1' => 'chrUn_gl000238',
		'GL000244.1' => 'chrUn_gl000244',
		'GL000248.1' => 'chrUn_gl000248',
		'GL000196.1' => 'chr8_gl000196_random',
		'GL000249.1' => 'chrUn_gl000249',
		'GL000246.1' => 'chrUn_gl000246',
		'GL000203.1' => 'chr17_gl000203_random',
		'GL000197.1' => 'chr8_gl000197_random',
		'GL000245.1' => 'chrUn_gl000245',
		'GL000247.1' => 'chrUn_gl000247',
		'GL000201.1' => 'chr9_gl000201_random',
		'GL000235.1' => 'chrUn_gl000235',
		'GL000239.1' => 'chrUn_gl000239',
		'GL000210.1' => 'chr21_gl000210_random',
		'GL000231.1' => 'chrUn_gl000231',
		'GL000229.1' => 'chrUn_gl000229',
		'GL000226.1' => 'chrUn_gl000226',
		'GL000207.1' => 'chr18_gl000207_random',
		'1' => 'chr1',
		'10' => 'chr10',
		'11' => 'chr11',
		'12' => 'chr12',
		'13' => 'chr13',
		'14' => 'chr14',
		'15' => 'chr15',
		'16' => 'chr16',
		'17' => 'chr17',
		'18' => 'chr18',
		'19' => 'chr19',
		'2' => 'chr2',
		'20' => 'chr20',
		'21' => 'chr21',
		'22' => 'chr22',
		'3' => 'chr3',
		'4' => 'chr4',
		'5' => 'chr5',
		'6' => 'chr6',
		'7' => 'chr7',
		'8' => 'chr8',
		'9' => 'chr9',
		'MT' => 'chrM',
		'X' => 'chrX',
		'Y' => 'chrY',
	};

	return ($chromNameHsh_ref);
}
sub logCalledCMDAndScript {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|114
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $resultScriptDir, $scriptAbsPath
#	output: 
#	toCall: &logCalledCMDAndScript($ARGVStr, $resultScriptDir, $scriptAbsPath);
#	calledInLine: 117
#....................................................................................................................................................#
	my ($ARGVStr, $resultScriptDir, $scriptAbsPath) = @_;


	my $cpScriptPath = "$resultScriptDir/script.ran.pl";
	my $calledCMDPath = "$resultScriptDir/called.cmd.txt";
	system "cp -f $scriptAbsPath $cpScriptPath";
	system "chmod 0444 $cpScriptPath"; #---[07/03/2014 18:02] make it read-only to make sure there'll be accodental change of parameters
	open CALLEDCMD, ">", $calledCMDPath;
	print CALLEDCMD join "", ($ARGVStr), "\n";
	close CALLEDCMD;
	
	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|276
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|73, 4_finishingTasks|114
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 76, 119
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->276
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->276
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->276
		print "=========================================================================\n\n";
	}
}
sub printOutputFileListAndReadme {
#....................................................................................................................................................#
#	subroutineCategory: output
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|114
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $outDir, $paramTag
#	output: 
#	toCall: &printOutputFileListAndReadme($ARGVStr, $paramTag, $outDir);
#	calledInLine: 118
#....................................................................................................................................................#
	my ($ARGVStr, $paramTag, $outDir) = @_;
	
	my $outputFileListPath = "$outDir/$paramTag/output.file.list.txt";
	open (OUTFILELIST, ">", $outputFileListPath);

	my %dirHsh = ();
	my %filelistLenCountHsh = ();
	push @{$filelistLenCountHsh{'dir'}}, length 'Directory';
	push @{$filelistLenCountHsh{'name'}}, length 'Name';
	push @{$filelistLenCountHsh{'description'}}, length 'Description';
	
	foreach my $outputFilePath (sort {$a cmp $b} keys %{$globalReadmeHsh_ref}) {
		my $fileDescription =  $globalReadmeHsh_ref->{$outputFilePath}{'description'};
		my $cleandOutputFilePath = $outputFilePath;
		$cleandOutputFilePath =~ s/\/+/\//g;
		
		my ($filePrefix, $fileDir, $fileSuffix) = fileparse($cleandOutputFilePath, qr/\.[^.]*/);
		$fileDir =~ s/^$outDir//;
		my $fileName = $filePrefix.$fileSuffix;
		$dirHsh{$fileDir}{$fileName} = $fileDescription;
		push @{$filelistLenCountHsh{'dir'}}, length $fileDir;
		push @{$filelistLenCountHsh{'name'}}, length $fileName;
		push @{$filelistLenCountHsh{'description'}}, length $fileDescription;
		
		open README, ">", "$outputFilePath.readme.txt";
		print README "=================\n";
		print README "File descriptions\n";
		print README "=================\n";
		print README "$fileDescription\n";
					
		if (exists $globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}) {
			my @colLenCountHsh = (length 'column');
			push @colLenCountHsh, length $_ foreach (@{$globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}});
			my $headerColLen = max(@colLenCountHsh)+2;
			print README "\n";
			print README "\n";
			print README "===================\n";
			print README "Column descriptions\n";
			print README "===================\n";
			print README "\n";
			printf README "%-".$headerColLen."s", 'column';
			print README "description\n";
			printf README "%-".$headerColLen."s", '------';
			print README "-----------\n";
			foreach my $header (@{$globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}}) {
				my $columnDescription = 'self-explanatory';
				$columnDescription = $globalReadmeHsh_ref->{$outputFilePath}{'header'}{$header} if exists $globalReadmeHsh_ref->{$outputFilePath}{'header'}{$header};
				printf README "%-".$headerColLen."s", $header;
				print README $columnDescription."\n";
			}
		}
		
		if (exists $globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}) {
			print README "\n";
			print README "\n";
			print README "=================\n";
			print README "Extra information\n";
			print README "=================\n";
			print README "\n";
			foreach my $title (sort keys %{$globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}}) {
				print README "$title\n";
				print README "-" foreach (1..length $title);
				print README "\n";
				print README "$_\n" foreach @{$globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}{$title}};
			}
		}
		
		print README "\n";
		print README "\n";
		print README "~" foreach (1..length "$fileName was created from running,");
		print README "\n";
		print README "$fileName was created from running,\n";
		print README "\n";
		print README "$ARGVStr\n";
		print README "\n";
		close README;
	}

	my $fileDir_colLen = max(@{$filelistLenCountHsh{'dir'}})+2;
	my $fileName_colLen = max(@{$filelistLenCountHsh{'name'}})+2;
	my $fileDescription_colLen = max(@{$filelistLenCountHsh{'description'}})+2;
	printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", 'directory', 'name', 'description');
	printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", '=========', '====', '===========');
	foreach my $fileDir (sort {$a cmp $b} keys %dirHsh) {
		foreach my $fileName (sort {$a cmp $b} keys %{$dirHsh{$fileDir}}) {
			my $fileDescription = $dirHsh{$fileDir}{$fileName};	
			printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", $fileDir, $fileName, $fileDescription);
		}
	}
	
	print OUTFILELIST "\n";
	print OUTFILELIST "\n";
	print OUTFILELIST "~" foreach (1..length "The above files were generated by running,");
	print OUTFILELIST "\n";
	print OUTFILELIST "The above files were generated by running,\n";
	print OUTFILELIST "\n";
	print OUTFILELIST "$ARGVStr\n";
	print OUTFILELIST "\n";

	close OUTFILELIST;

	return ();
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|73
#	secondaryAppearInSection: >none
#	input: none
#	output: $bamPath, $maxLen, $minLen, $outDir, $outputPrefix
#	toCall: my ($bamPath, $minLen, $maxLen, $outputPrefix, $outDir) = &readParameters();
#	calledInLine: 77
#....................................................................................................................................................#
	
	my ($bamPath, $minLen, $maxLen, $exclude_flag, $longShortFormat, $min_MAPQ, $outputPrefix, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/output/";
	$longShortFormat = 'short';
	
	GetOptions 	(
		"bamPath=s"  => \$bamPath,
		"minLen:i"  => \$minLen,
		"maxLen:i"  => \$maxLen,
		"min_MAPQ:i"  => \$min_MAPQ,
		"exclude_flag:s"  => \$exclude_flag,
		"longShortFormat:s"  => \$longShortFormat,
		"outputPrefix=s"  => \$outputPrefix,
		"outDir:s"  => \$outDir
	)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($bamPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 777 $outDir/";
	
	return($bamPath, $minLen, $maxLen, $exclude_flag, $longShortFormat, $min_MAPQ, $outputPrefix, $outDir);

}
sub reportAndLogStatus {
#....................................................................................................................................................#
#	subroutineCategory: log
#	dependOnSub: currentTime|276
#	appearInSub: checkBedtoolsVersion|158, checkSamtoolsVersion|178, convertBamToCTSS|198
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|73, 3_convert|105
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportAndLogStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 172, 192, 213
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->276
	print TMPLOG "[".&currentTime()."] ".$message.$lineEnd if $lineEnd ne "\r";#->276
	
	return ();
}
sub timeStamp {
#....................................................................................................................................................#
#	subroutineCategory: time, general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: none
#	output: $curntTimeStamp
#	toCall: my ($curntTimeStamp) = &timeStamp();
#	calledInLine: 60
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $curntTimeStamp = sprintf "%04d.%02d.%02d.%02d.%02d.%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec;	

	return ($curntTimeStamp);
}

exit;


















































