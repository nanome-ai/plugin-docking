# Nanome - Docking

A Nanome Plugin to interface with a variety of docking softwares to dock ligands to a receptor.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

To use Autodock4 on Windows, Autodock4 must be installed on the computer and in the PATH variable.

## Usage

To run Docking Smina in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address> [args]
```

To run Docking Autodock4 on Windows:

```
$ python3 -m pip install -r requirements.txt
$ python3 run.py -a <plugin_server_address> autodock4 [args]
```

Note: requires Autodock4 to be installed on the computer and in the PATH variable.

---

In Nanome:

- Activate Plugin, Docking window will automatically open
- Select a receptor
- Click on "Ligand", and select ligands to dock
- If using Smina, click on "Site", and select which molecule should be used to define the docking site
- Choose exhaustiveness, number of results and size of the box to generate around the site molecule (Smina only)
- If visual scoring is turned on, atom size and labels will indicate each atom's contribution to the ligand's score
- Click on Run

## Development

To run Docking with autoreload:

```
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> smina [args]
```

## License

MIT
