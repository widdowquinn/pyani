FROM python:3
MAINTAINER Leighton Pritchard <leighton.pritchard@hutton.ac.uk>

RUN apt-get update && apt-get install -y \
			       ncbi-blast+ \
			       mummer \
			       && \
    pip3 install --upgrade pip && \
    pip3 install pyani

WORKDIR /host_dir

ENTRYPOINT ["average_nucleotide_identity.py"]