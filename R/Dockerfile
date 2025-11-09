FROM rocker/r-ver:4.5.1

# make renv install destination persistent at runtime
ENV RENV_PATHS_LIBRARY=/usr/local/lib/R/site-library

WORKDIR /project

COPY renv/ /project/renv/
COPY renv.lock .Rprofile master-thesis.Rproj /project/

RUN apt-get update && apt-get install -y libicu-dev pandoc \
  && R -e "install.packages('renv', repos='https://cloud.r-project.org')" \
  && mkdir -p /usr/local/lib/R/site-library \
  && R -e "source('/project/renv/activate.R'); renv::restore(prompt = FALSE, clean = TRUE, rebuild = TRUE)"

COPY . /project

CMD ["R"]
