# CAGE_HCC
Analysis of CAGE data for liver cancer samples

How to process the CAGE data from BAM to output analysis

# STEP 1: BAM TO CTSS
It generated one CTSS file for one BAM file (one sample). The script is BamtoCTSS.pl. 
Original path: /home/hon-chun/resources/perlScript/FANTOM/bamToCTSS/v0.6/bamToCTSS_v0.6.pl

Wrapper bash script to submit bamToCTSS_v0.6.pl as a job on the cluster for one or several samples:

BamToCTSS.sh
  usage: 
		./BamToCTSS.sh input.list output_path
  arguments: 
		input.list: list of paths pointing to BAM files
  	output_path: path to output directory where to save the CTSS files


#work in progress...
#STEP 2: mergeCTSSfile

We gotta pool the CTSS before calling paraclu clusters.
        For that purpose I use this script:

        /home/hon-chun/resources/perlScript/FANTOM/CTSSBedPooler/v0.2/CTSSBedPooler_v0.2.pl

        For the command of running CTSSBedPooler_v0.2.pl, check here:
        /home/hon-chun/resources/perlScript/FANTOM/CTSSBedPooler/v0.2/CTSSBedPooler_v0.2.cmd.log.txt

        It is a simple wrapper of the bedtools merge command:

        gzip -dc /path/to/dir/*ctss.bed.gz | bedtools merge -s -d -1 -i stdin -c 4,5,6 -o sum,sum,distinct | gzip -c >/path/to/pooled.ctss.bed.gz

As the paraclu manual describes the merged input files to paraclu need to be sorted by strand before position, i.e all + reads for chr1 appear together.


STEP 3: PARACLU
http://cbrc3.cbrc.jp/~martin/paraclu/
Define CAGE clusters: paraclu. Parameters: minValue per cluster = 30. Subset the output of paraclu: paraclu-cut. Default parameters: length of cluster <200; density/baseline density <2

Other: To call enhancers from CAGE data: https://github.com/anderssonrobin/enhancers
