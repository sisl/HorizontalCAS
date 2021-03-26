FROM julia:1.1.1
WORKDIR /code
RUN apt-get update
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
    Pkg.PackageSpec(;name="Interact", version="v0.10.3"), \
    Pkg.PackageSpec(;name="PGFPlots", version="v3.2.1"), \
    Pkg.PackageSpec(;name="Colors", version="v0.12.6"), \
    Pkg.PackageSpec(;name="ColorBrewer", version="v0.4")])'
