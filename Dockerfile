FROM julia:1.1.1
WORKDIR /code
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y --fix-missing --no-install-recommends \
    debconf build-essential autoconf \
    curl dirmngr apt-transport-https lsb-release ca-certificates \
    python3 python3-setuptools python3-pip
# install latex dependencies
# if timezone is not setup texlive-latex-base can not be installed
RUN ln -snf /usr/share/zoneinfo/Etc/UTC /etc/localtime \
    && echo "Etc/UTC" > /etc/timezone
RUN apt-get install -y --fix-missing --no-install-recommends \
    pdf2svg texlive-pictures texlive-latex-extra texlive-luatex
# install nodejs
RUN curl -sL https://deb.nodesource.com/setup_14.x | bash -
RUN apt-get install -y --fix-missing --no-install-recommends nodejs
RUN python3 -m pip -q install pip --upgrade
RUN python3 -m pip install jupyterlab==2.3.1
# RUN useradd -ms /bin/bash jupyter
# USER jupyter
RUN julia -e 'using Pkg; Pkg.add.([ \
    Pkg.PackageSpec(;name="Printf"), \
    Pkg.PackageSpec(;name="POMDPs", version="v0.7.0"), \
    Pkg.PackageSpec(;name="POMDPModelTools", version="v0.1.2"), \
    Pkg.PackageSpec(;name="LocalFunctionApproximation", version="v1.1.0"), \
    Pkg.PackageSpec(;name="GridInterpolations", version="v1.1.1"), \
    Pkg.PackageSpec(;name="Distributed"), \
    Pkg.PackageSpec(;name="SharedArrays"), \
    Pkg.PackageSpec(;name="StaticArrays", version="v0.12.5"), \
    Pkg.PackageSpec(;name="HDF5", version="v0.12.5"), \
    Pkg.PackageSpec(;name="Revise"), \
    Pkg.PackageSpec(;name="Interact", version="v0.10.3"), \
    Pkg.PackageSpec(;name="PGFPlots", version="v3.2.1"), \
    Pkg.PackageSpec(;name="Colors", version="v0.12.6"), \
    Pkg.PackageSpec(;name="ColorBrewer", version="v0.4"), \
    Pkg.PackageSpec(;name="IJulia"), \
    Pkg.PackageSpec(;name="WebIO")])'
ENV JUPYTER_ENABLE_LAB=yes
RUN julia -e 'using WebIO; WebIO.install_jupyter_labextension()'
EXPOSE 8888/tcp
