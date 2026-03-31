FROM julia:1.10-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends curl wget git build-essential ca-certificates libopenblas-dev liblapack-dev libhdf5-dev hdf5-tools python3 python3-pip python3-venv libgl1 libglib2.0-0 libxrender1 libxext6 libx11-6 libxt6 libxft2 htop tree ripgrep jq && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 -s /bin/bash astro
USER astro
WORKDIR /home/astro/project

ENV JULIA_PROJECT=/home/astro/project
ENV JULIA_DEPOT_PATH=/home/astro/.julia
ENV JULIA_NUM_THREADS=auto
ENV JULIA_REVISE_INCLUDE=1
ENV GKSwstype=100

COPY --chown=astro:astro julia/Project.toml ./Project.toml

RUN julia --project=. -e 'import Pkg; Pkg.add(["SatelliteToolboxBase","SatelliteToolboxPropagators","SatelliteToolboxTle","SatelliteToolboxTransformations","OrdinaryDiffEq","StaticArrays","Plots","CSV","DataFrames","Revise","BenchmarkTools","ProgressMeter"]); println("OK")'

COPY --chown=astro:astro julia/src ./src
COPY --chown=astro:astro julia/scripts ./scripts

EXPOSE 1234 8888

CMD ["julia", "--project=.", "-i", "-e", "using Revise; println(\"Astrodynamics ready!\")"]
