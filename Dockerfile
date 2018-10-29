FROM marsluo/ubuntu_r3.5.1

# link 'Rscript'
RUN sudo ln -s /opt/R/3.5.1/bin/Rscript /usr/bin/Rscript

# install R packages: shiny, shinyjs, shinythemes..
RUN sudo R -e "install.packages('shiny', repos='https://cran.rstudio.com/')"
RUN sudo R -e "install.packages('shinyjs', repos='https://cran.rstudio.com/')"
RUN sudo R -e "install.packages('shinythemes', repos='https://cran.rstudio.com/')"


RUN sudo apt-get update
RUN sudo apt-get -y upgrade

# install gdebi-core
RUN sudo apt-get -y install gdebi-core
RUN sudo wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb
RUN sudo gdebi -n shiny-server-1.5.9.923-amd64.deb

# install openssl
RUN sudo apt-get install -y libssl-dev

# install libxml2-dev
RUN sudo apt-get install -y libxml2-dev


# install required R packages for shinyapps of 2nd analysis
# RUN sudo R -e "r <- getOption(\"repos\");r[\"CRAN\"] <- \"https://cran.rstudio.com/\";r[\"BioCsoft\"] <- \"https://bioconductor.org/packages/3.8/bioc\";r[\"BioCann\"] <- \"https://bioconductor.org/packages/3.8/data/annotation\";r[\"BioCexp\"] <- \"https://bioconductor.org/packages/3.8/data/experiment\";r[\"BioCworflows\"] <- \"https://bioconductor.org/packages/3.8/workflows\";options(repos = r)"

RUN sudo R -e "install.packages('XML',repos='https://cran.rstudio.com/')"

RUN sudo R -e "install.packages('Rtsne',repos='https://cran.rstudio.com/')"

RUN sudo R -e "install.packages('pheatmap',repos='https://cran.rstudio.com/')"

# Copy Shiny application into Docker container
COPY UI.R /srv/shiny-server/UI.R
COPY Server.R /srv/shiny-server/Server.R

# Copy Shiny configuration into Docker container
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /srv/shiny-server/shiny-server.sh

# Copy data into Docker container
COPY GRCh38/barcodes.tsv /srv/shiny-server/GRCh38/barcodes.tsv
COPY GRCh38/genes.tsv /srv/shiny-server/GRCh38/genes.tsv
COPY GRCh38/matrix.mtx /srv/shiny-server/GRCh38/matrix.mtx

# Copy CSS styling and ggplot theming into Docker container


EXPOSE 3838
EXPOSE 80

# link 'shiny-server.sh'
RUN sudo ln -s /srv/shiny-server/shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]