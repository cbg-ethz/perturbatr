FROM r-base

RUN mkdir /analysis/
RUN mkdir /analysis/data/
RUN mkdir /analysis/plots/

RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('dtplyr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('tidyr')"
RUN Rscript -e "install.packages('lme4')"
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('ggthemr')"
RUN Rscript -e "install.packages('hrbrthemes')"
RUN Rscript -e "install.packages('viridis')"
RUN Rscript -e "install.packages('cowplot')"
RUN Rscript -e "install.packages('formula.tools')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('xtable')"

RUN Rscript -e "devtools::install_github('cbg-ethz/perturbatr@v0.99.0')"

COPY run.sh /analysis/
COPY *.R /analysis/
COPY data/rnai_screen_raw.rds /analysis/

CMD ["/analysis/run.sh"]
