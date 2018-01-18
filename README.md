# Clonal Growth Simulator (ABM)

- [Installation](#installation)
- [Usage](#usage)
  - [Run](#run)
  - [Configuration files](#configuration-files)
- [Docker](#docker)
- [Acknowledgements](#acknowledgements)


## Installation

Dependencies:
- A C++ compiler that supports C++11
- Boost libraries: system, filesystem, iostreams, random
- CMake 2.6 or newer
- Gzip (optional)

On Ubuntu 16.04 all dependencies can be installed with the following command:

``` install cmake libboost-filesystem-dev libboost-iostreams-dev libboost-random-dev gcc g++ python gzip```

Build:
```bash
mkdir build
cd build
cmake ../
make
```


## Usage

### Run
To run a simulation you need to call the binary with an xml file:

```bash
mkdir out
./bin/simulator examples/iteratedgrowth_simple.xml
```

### Configuration files
The simulation settings are stored in a single xml file containing 2 nodes: `Simulation` and `Experiment`.
An example of such a file can be found in the `examples/` folder.


#### Model
The model node specifies the growth model parameters:

| Element        | Description                                 |
|----------------|---------------------------------------------|
| DivisionRate   | Number of times cells divide per day        |
| DeathRate      | Number of times cells die per day        |
| MutationSD     | SD of the normal distribution used for mutations |
| DivisionRateSD | SD of the normal distribution used to set the initial division rates |
| MinDivisionRate | minimum division rate |
| MaxDivisionRate | maximum division rate |

#### Experiment

The experiment node specifies the kind of experiments and the experiment specific settings. The experiment is specified via the `type` attribute, and can be:

- **ConstantGrowth**: cells grow for a given time
- **IteratedGrowth**: cells grow until a target number is reached and then a subset is passed to the next generation.

The following parameters can be used for both experiments:

| Element               | Description                                                                    |
|-----------------------|--------------------------------------------------------------------------------|
| Name                  | Simulation name, used for filenames of the output                              |
| Outpath               | Path to store output to                                                        |
| SaveXML               | Export all settings (including defaults) in an xml file and save it to Outpath |
| gzip                  | Gzip all output files                                                          |
| UseSimDir             | Create a folder in Outpath that contains all generated output                  |
| TrackDivTime          | Save all division time to file |
| TrackCells            | Save cells |
| InitialPopulationSize | Initial number of cells                                                        |
| InitUniform           | Number of clones over which the initial cells are uniformly distributed        |
| InitFile              | File with initial clone distribution                                           |
| InitSeedDistr         | Seed used for setting up the master population                                 |
| InitSeed              | Seed used for taking a sample form the master population                   |
| Seed                  | Seed used for stochastic growth and passage                                    |


##### Parameters for ConstantGrowth

| Element               | Description                                                                    |
|-----------------------|--------------------------------------------------------------------------------|
| SimulationTime        | Simulation time in days                                                        |
| SaveInterval          | Time interval for saving                                                |

##### Parameters for IteratedGrowth

| Element                  | Description                                                                    |
|--------------------------|--------------------------------------------------------------------------------|
| CriticalPopulationSize   | Population size after which a growth step is stopped                           |
| NumberOfCellsToKeep      | Number of cells to keep after passage                                          |
| NumberOfPassages         | Number of growth and passage cycles                                            |
| MaxPassTime              | Maximum time a growth step may take                                            |
| StopSimIfNPassNotReached | Stop the simulation if MaxPassTime is reached                                  |


## Docker
Alternatively, you can use Docker to build the software:

```docker build -t abm .```

Then, you can run the software with:

```
cd examples/
mkdir out
docker run -u $(id -u):$(id -g) -v $PWD:/data/ abm /data/it_abm_docker.xml
```

Note that this command mounts the current directory as `/data/` inside the container. Therefore, `OutPath` in the xml must start with `/data/` as well.

## Acknowledgements

We thank the developers of the TinyXML-2 library for providing this library (https://github.com/leethomason/tinyxml2).
