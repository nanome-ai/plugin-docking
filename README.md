# Nanome - Docking

A Nanome Plugin to interface with a variety of docking softwares to dock ligands to a receptor.

Supported Docking algorithms:
- Smina
- Autodock4
- Rhodium (Coming soon)

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

## Usage
The docking algorithm used is determined by the `--algorithm` flag on startup. It flag not provided, it defaults to `smina`

To run Smina Docking in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh [args]
```

To run Docking with Autodock4:
```
$ cd docker
$ ./build.sh --build-arg ALGORITHM=autodock4
$ ./deploy.sh --algorithm autodock4 -a [args]
```
---

In Nanome:

- Activate Plugin, Docking window will automatically open
- Select a receptor
- Click on "Ligand", and select ligands to dock
- If using Smina, click on "Site", and select which molecule should be used to define the docking site
- Choose number of poses to return for each ligand
- Choose size of the box to generate around the site molecule (Smina only)
- If visual scoring is turned on, atom size and labels will indicate each atom's contribution to the ligand's score
- Click Run

## Development

To run Docking with autoreload:
```
$ python3 -m pip install -r requirements.txt
$ python3 run.py <algorithm> -r [args]
```
algorithm can be (smina | autodock4)

Note for autodock4: The adfr-suite conda environment must be set up to run code requiring Python 2.7
```sh
$ conda env create --file adfr-suite.yml
```

## License

MIT
