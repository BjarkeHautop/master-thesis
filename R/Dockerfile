# Use R 4.5.1 base image
FROM rocker/r-ver:4.5.1

# Set working directory inside container
WORKDIR /project

COPY renv/ /project/renv/
COPY renv.lock .Rprofile master-thesis.Rproj /project/

# Install renv and restore packages non-interactively
RUN apt-get update && apt-get install -y libicu-dev pandoc \
    && R -e "install.packages('renv', repos='https://cloud.r-project.org')" \
    && R -e "source('renv/activate.R'); renv::restore(prompt = FALSE, clean = TRUE)"

# Copy the rest of the project (source code, data, etc.)
COPY . /project

# Default: launch R (can override with docker run)
CMD ["R"]
