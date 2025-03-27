FROM python:3.10-slim-buster
LABEL maintainer=shuzhao

RUN apt install apt-transport-https dirmngr gnupg ca-certificates; \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF; \
    echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list; \
    apt update; \
    apt install -y --no-install-recommends wget unzip mono-devel; \
    pip install --upgrade numpy; \
    pip install --upgrade asari-metabolomics; \
    mkdir -p /usr/local/thermo/; \
    wget -O /usr/local/thermo/ThermoRawFileParser.zip "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.0/ThermoRawFileParser.zip"; \
    cd /usr/local/thermo/; \
    unzip  ThermoRawFileParser.zip; \
    rm -rf ThermoRawFileParser.zip; \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["bin/bash"]

