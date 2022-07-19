# syntax=docker/dockerfile:1

# use 
FROM ubuntu:jammy

# install some packages
RUN apt update
RUN apt-get install -y python3.9 python3-pip python3-setuptools python3-dev python3-parasail
RUN apt-get install -y ncbi-blast+ fasttree
RUN apt-get install -y flex bison libgmp3-dev cmake mpich # for compiling generax
RUN apt-get -y install default-jre # for astral
RUN apt-get install -y nano wget git unzip

# install required python modules
RUN pip3 install numpy scipy Bio

# create the exec directory to hold muscle, fasttree, and generax
RUN mkdir /exec

# install muscle
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
RUN chmod ugo+x muscle5.1.linux_intel64
RUN mv muscle5.1.linux_intel64 /exec/muscle

# install generax
RUN git clone --recursive https://github.com/BenoitMorel/GeneRax
#WORKDIR "GenRax"
RUN cd /GeneRax && sh install.sh

# install ASTRAL
RUN wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
RUN unzip Astral.5.7.8.zip
RUN rm Astral.5.7.8.zip

# set up symbolic links to /usr/local/bin
RUN ln -s /exec/muscle /usr/local/bin/muscle
RUN ln -s /usr/bin/fasttree /usr/local/bin/fasttree
RUN ln -s /GeneRax/build/bin/generax /usr/local/bin/generax
RUN ln -s /Astral/astral.5.7.8.jar /usr/local/bin/astral

# clone xenoGI
RUN git clone -b devel https://github.com/ecbush/xenoGI.git


