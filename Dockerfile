FROM ubuntu:18.04

MAINTAINER Jacob Pfeil, jpfeil@ucsc.edu

# Update and install required software
RUN apt-get update --fix-missing

RUN apt-get install -y wget default-jre

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda 

ENV PATH=/root/miniconda/bin:$PATH

RUN conda update -y conda && conda install -y -c bioconda samtools bowtie2 bamtools

WORKDIR /opt

# Add Trimmomatic
COPY Trimmomatic-0.39 /opt/trim

# Add wrapper scripts
COPY pipeline /opt/pipeline

# Data processing occurs at /data
WORKDIR /data

ENTRYPOINT ["python", "/opt/pipeline/run.py"]
CMD ["-h"]
