#!/bin/sh

#Description
#*-*-*-*-*-*-
# This script submit a job to run bamToCTSS_v0.6.pl
# which convert a BAM file from CAGE to a CTSS FILE
# It takes two parameters: (1) a text file with the path of BAM file(s), (2) the output directory for the CTSS files

#------------
#

#Usage
#*-*-*-
#bash ./BamToCTSS.sh input.list output_path
#-----

#Variables
#*-*-*-
INPUT=${1}
OUTDIR=${2}
SCRIPT=/home/pauline/scripts/CAGE_processing/bamToCTSS_v0.6.pl
TEMP=/home/pauline/scripts/tmp
#------


#code
#*-*-
mkdir -p ${OUTDIR}

#create one bash script per sample
cat ${INPUT} | while read line
do
  DATAID=$(basename ${line} .bam)
  echo " processing $DATAID"

  printf "#!/bin/sh\nperl ${SCRIPT} --bamPath=${line} --minLen=0 --maxLen=100 --min_MAPQ=10 --longShortFormat=short --exclude_flag=512,256,1024,2048 --outputPrefix=$DATAID --outDir=${OUTDIR}\ngunzip -k ${OUTDIR}/${DATAID}/bed/${DATAID}.short.ctss.bed.gz" > ${TEMP}/${DATAID}.BamToCTSS.sh

	qsub -o ${TEMP}/${DATAID}.BamToCTSS.sh.o -e ${TEMP}/${DATAID}.BamToCTSS.sh.e ${TEMP}/${DATAID}.BamToCTSS.sh


#logs
cp ${TEMP}/${DATAID}.BamToCTSS.sh ${OUTDIR}/${DATAID}/bed/readme_${DATAID}.BamToCTSS.sh.txt
echo -e "LOGS \nINPUT=${1}\nOUTDIR=${2}\nSCRIPT=${SCRIPT}\nTEMP=${TEMP}" >> ${OUTDIR}/${DATAID}/bed/readme_${DATAID}.BamToCTSS.sh.txt


done


#------

#logs
#*-*-
cp /home/pauline/scripts/CAGE_processing/BamToCTSS.sh ${OUTDIR}/readme_BamToCTSS.sh.txt
echo -e "LOGS \nINPUT=${1}\nOUTDIR=${2}\nSCRIPT=${SCRIPT}\nTEMP=${TEMP}" >> ${OUTDIR}/readme_BamToCTSS.sh.txt
cp ${OUTDIR}/readme_BamToCTSS.sh.txt /home/pauline/scripts/by_date_copy_script_run/$(date -d "today" +"%Y%m%d%H%M")_BamToCTSS.sh
