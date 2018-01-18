FROM ubuntu:16.04
MAINTAINER Margriet Palm

RUN apt-get update -y && \
    apt-get install -y cmake libboost-filesystem-dev \
                       libboost-iostreams-dev libboost-random-dev \
                       gcc g++ python gzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . /src/ClonalGrowthSimulator_ABM

# build ClonalGrowthSimulator_ABM
RUN mkdir -p ClonalGrowthSimulator_ABM/build/
RUN cd ClonalGrowthSimulator_ABM/build/ && cmake ../ && make

WORKDIR /app
RUN mv /src/ClonalGrowthSimulator_ABM/bin/simulator /app/simulator
RUN rm -rf /src

ENTRYPOINT ["/app/simulator"]
