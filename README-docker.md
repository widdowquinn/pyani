# README-docker.md

This document describes using building/testing a Docker container for `pyani`

## Set up the registry

Start the registry container.

```
$ docker run -d -p 5000:5000 --restart=always --name registry registry:2
```

## Build the image

```
$ docker build -t leightonpritchard/pyani .
```

Note on passing arguments: https://stackoverflow.com/questions/32727594/how-to-pass-arguments-to-shell-script-through-docker-run

Note on Python image: https://store.docker.com/images/python