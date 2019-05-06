# Generated by Dockter 0.14.4 at 2019-05-06T17:35:44.671Z
# To stop Dockter generating this file and start editing it yourself,
# rename it to "Dockerfile".

# This tells Docker which base image to use.
FROM ubuntu:18.04

# This section installs system packages needed to add extra system repositories.
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      apt-transport-https \
      ca-certificates \
      curl \
      software-properties-common

# This section adds system repositories required to install extra system packages.
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-add-repository "deb https://mran.microsoft.com/snapshot/2019-05-05/bin/linux/ubuntu bionic-cran35/"

# This section sets environment variables within the image.
ENV TZ="Etc/UTC"

# This section installs system packages required for your project
# If you need extra system packages add them here.
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      git-core \
      libcurl4-openssl-dev \
      libssh2-1-dev \
      libssl-dev \
      libxml2-dev \
      make \
      pandoc pandoc-citeproc \
      r-base \
      zlib1g-dev \
 && apt-get autoremove -y \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# It's good practice to run Docker images as a non-root user.
# This section creates a new user and its home directory as the default working directory.
RUN useradd --create-home --uid 1001 -s /bin/bash dockteruser
WORKDIR /home/dockteruser

# This is a special comment to tell Dockter to manage the build from here on
# dockter

# This section copies package requirement files into the image
COPY .DESCRIPTION DESCRIPTION

# This section runs commands to install the packages specified in the requirement file/s
RUN bash -c "Rscript <(curl -sL https://unpkg.com/@stencila/dockter/src/install.R)"

# This sets the default user when the container is run
USER dockteruser