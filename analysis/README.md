# Data analysis

The folder contains files needed to recreate the data analysis for the presented paper.

## Usage

In order to guarantee best reproducibility of our results, we prepared a Docker image. Make sure to have Docker installed and the daemon running
using:

```bash
  docker
```

If this does not succeed, install Docker from [here](https://www.docker.com/community-edition).

Afterwards load the image with:

```bash
  docker load < docker_image.tar
```

The image can then be used with:

```bash
  docker run analysis
```

This opens a terminal session from which you can call:

```bash
  /analysis/run.sh
```

in order to redo the analysis.
