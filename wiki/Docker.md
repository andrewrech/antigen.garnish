# Docker instructions.

`antigen.garnish` can be non-interactively run using a public docker image in a local container. A large wrapper script takes vcf input, HLA table or JSON input, and generates all possible output types. Customizability is limited when running non-interactively. Only direct VCF input is supported in this mode. An additional docker image is available to annotate vcfs using [SnpEff](https://sourceforge.net/projects/snpeff/).

## Install Docker
Install [Docker](https://www.docker.com/get-started) and start the Docker daemon. Set appropriate memory limits and core usage for your machine in the Docker app and then open up a terminal (on Mac) or the Linux Virtual Machine command line in the Docker Windows GUI.

## Running antigen.garnish with Docker

### pull the antigen.garnish docker image
```sh

docker pull leeprichman/antigen_garnish

```

### snpEff: Annotate your vcfs with [SnpEff](https://sourceforge.net/projects/snpeff/)
```sh
# start up our container
cID=$(docker run -it -d leeprichman/snpeff /bin/bash)

# set up and move our input file from working directory (change names as appropriate here)
# one single sample vcf file at a time
VCF="myvcf.vcf"

# copy our files on
docker cp $VCF $cID:/$VCF

# run the script to hg19 annotate (substitute others here. See SnpEff databases)
# GRCm38.86 and hg38 are preloaded on this docker image, snpEff will download others
docker exec $cID snpeff.sh hg19

# output file name will be with _se inserted
VCFO=$(echo $VCF | sed 's/\.vcf/_se\.vcf/')
# copy output back to the working directory into output folder
docker cp $cID:/$VCFO .

# clean up the container for next sample
docker stop $cID
docker rm  $cID
```

### antigen.garnish: Start container, move files on, execute wrapper script, and recover output
```sh
# start up our container
cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)

# set up and move our input files from working directory (change names as appropriate here)
# one single sample, snpEff annotated vcf file at a time
# for tumor allelic fraction to be properly recognized the vcf name must be the same as the vcf sample column
# for many vcfs, the column names are "TUMOR", "NORMAL"
# TUMOR.vcf, TUMOR_se.vcf, TUMOR.ann.vcf, and TUMOR.vcf.gz all work
# MHC may be a .JSON file as from xHLA output or a 1 row csv or tsv named "*mhc.txt" formatted like:
#       sample_id     MHC
#       mysample.vcf  "H-2-Kb H-2-Db H-2-IAb"
# MHC must be a quoted space-separated string of MHC alleles in antigen.garnish format (see ?list_MHC)
# optional RNA transcript level count matrix to be matched against (minimum tpm = 1), named "*counts*"
# must contain a first column of Ensembl Transcript IDs, and a column named "tpm".
VCFO="myvcf_se.vcf"
MHC="mhc.json"
# optional
RNA="rna_counts.txt"

# copy our files on
docker cp $VCFO $cID:/$VCFO
docker cp $MHC $cID:/$MHC
# optional
docker cp $RNA $cID:/$RNA

# run the big  wrapper script
docker exec $cID run_antigen.garnish.R

# copy output back to the working directory into output folder
docker cp $cID:/ag_docker_output/ .

# clean up the container for next sample
docker stop $cID
docker rm  $cID
```
