FROM ubuntu

WORKDIR /single-cell

ENV DEBIAN_FRONTEND=noninteractive

COPY . .

RUN apt update \
 && apt -y upgrade \
 && apt install -y python3-pip python3-dev openjdk-8-jre \
    libvcflib-tools libvcflib-dev libcurl4-openssl-dev libffi-dev  \
    build-essential wget curl vim git software-properties-common samtools \
 && cp /usr/bin/python3 /usr/bin/python \
 && python -m pip install -r requirements.txt \ 
 && chmod 777 ./software/parallel \
 && chmod 777 ./software/freebayes \
 && apt-get purge -y r-base* r-recommended r-cran-* \
 && apt -y autoremove \ 
 && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
 && apt update \
 && apt install -y r-base r-base-core r-recommended r-base-dev \
 && apt install -y libcurl4-gnutls-dev libxml2-dev libssl-dev \
 && Rscript requirements.R \
 && chmod 777 software/* \
 && curl -L -o ./software/structure_kernel_source.tar.gz \
	https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_kernel_source.tar.gz \
 && cd ./software \
 && tar zxvf structure_kernel_source.tar.gz \
 && cd structure_kernel_src/ \
 && make \
 && cp structure /usr/local/bin/ \
 && rm -rf /var/chache/apk/*

ENTRYPOINT ["snakemake"]

