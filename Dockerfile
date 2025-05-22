FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install only what you need
RUN Rscript -e "BiocManager::install('PharmacoGx', ask = FALSE, update = FALSE)"

# Optional: reduce image size by stripping debug symbols from native libs
RUN strip /usr/local/lib/R/site-library/*/libs/*.so || true

# Set working directory
WORKDIR /app

# Copy your scripts (optional)
COPY . .

# Default command — can override with `docker run ... Rscript your_script.R`
CMD ["Rscript"]