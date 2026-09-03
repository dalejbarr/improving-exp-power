FROM rocker/tidyverse:4.5.3

RUN wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
RUN tar xvzf gsl-latest.tar.gz
RUN cd gsl-2.8 && ./configure && make && make install

RUN Rscript -e "install.packages(c('remotes', 'devtools', 'binom'), repos='https://cran.rstudio.com/')"
# RUN Rscript -e "remotes::install_github('dalejbarr/explan')"

COPY psrsim/ /var/psrsim/
COPY manuscript/dotemacs /root/.emacs

RUN Rscript -e "devtools::install('/var/psrsim')"

RUN sudo apt update && apt-get install -y emacs-nox texlive-latex-base texlive-latex-recommended texlive-latex-extra

COPY scripts/ /var/scripts/
COPY results/ /var/results/
COPY csl/ /var/csl/

RUN chmod u+x /var/scripts/emacs-pkg-install.sh
RUN /var/scripts/emacs-pkg-install.sh ess
RUN /var/scripts/emacs-pkg-install.sh citeproc

COPY docker_int_do-all-simulations.sh /usr/local/bin/do-all-simulations.sh
RUN chmod a+x /usr/local/bin/do-all-simulations.sh
COPY docker_int_choose-num-basis-full.sh /usr/local/bin/do-choose-num-basis-full.sh
RUN chmod a+x /usr/local/bin/do-choose-num-basis-full.sh
COPY docker_int_make-manuscript.sh /usr/local/bin/make-manuscript
RUN chmod a+x /usr/local/bin/make-manuscript
COPY docker_int_make-supplementary.sh /usr/local/bin/make-supplementary
RUN chmod a+x /usr/local/bin/make-supplementary

WORKDIR "/var"

RUN mkdir /output
COPY manuscript/ /var/manuscript/
COPY supplementary /var/supplementary/
