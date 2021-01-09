FROM ubuntu:19.04
MAINTAINER Xiaoli Dong <xiaolid@gmail.com>	
LABEL version="1.0.4"
	
WORKDIR /NGStools/

#Install compiler and perl stuff
RUN apt-get update && apt-get install -y \
    #apt-utils \
    autoconf \
    #build-essential \
    cpanminus \
    gcc-multilib \
    git \
    make \
    openjdk-8-jdk \
    perl \
    python \
    sqlite3 \
    tar \
    unzip \
    wget 
    

# Install libraries that BioPerl dependencies depend on
RUN apt-get update && apt-get install -y \
    expat \
    graphviz \
    libdb-dev \
    libgdbm-dev \
    libexpat1 \
    libexpat-dev \
    libssl-dev \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev

#install perl modules
RUN cpanm Bio::Perl \
    DBI \
    Archive::Extract \
    DBD::SQLite \
    File::Copy::Recursive \
    Bio::DB::EUtilities \
    LWP::Protocol::https && \
    git clone https://git.code.sf.net/p/swissknife/git swissknife-git && \
    cd swissknife-git && \
    perl Makefile.PL && \
    make install && \
    cd /NGStools

#aragorn
RUN git clone https://github.com/TheSEED/aragorn.git && \
    cd aragorn && \
    gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.36.c && \
    cd /NGStools

#hmmer rRNAFinder need it
RUN git clone https://github.com/EddyRivasLab/hmmer && \
    cd hmmer && \
    git clone https://github.com/EddyRivasLab/easel && \
    autoconf && \
    ./configure && \
    make  && \
    cd /NGStools

#blast for classifying rRNA sequences
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.9.0+-x64-linux.tar.gz && \
    cd /NGStools

#prodigal
RUN git clone https://github.com/hyattpd/Prodigal.git && \
    cd Prodigal && \
    make && \
    cd /NGStools

#minced
RUN git clone https://github.com/ctSkennerton/minced.git && \
    cd minced && \
    make && \
    cd /NGStools

#diamond
RUN mkdir diamond && \
   cd diamond && \
    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && \
    rm diamond-linux64.tar.gz diamond_manual.pdf && \
    cd /NGStools

#MinPath
RUN wget http://ebg.ucalgary.ca/metaerg/minpath1.4.tar.gz && \
    tar -xzf minpath1.4.tar.gz && \
    rm minpath1.4.tar.gz && \
    cd /NGStools

#metaerg
RUN git clone https://github.com/xiaoli-dong/metaerg.git

# Clean
RUN apt-get remove -y autoconf \
    cpanminus \
    gcc-multilib \
    git \ 
    make && \
    apt-get autoclean -y 

ENV MinPath /NGStools/MinPath
ENV PATH="/NGStools/aragorn:/NGStools/minced:/NGStools/Prodigal:/NGStools/ncbi-blast-2.9.0+/bin:/NGStools/diamond:/NGStools/hmmer/src:/NGStools/MinPath:/NGStools/metaerg/bin:${PATH}"

WORKDIR /NGStools/