# Nanome - Docking

This Nanome plugin interfaces with a variety of docking softwares to dock ligands to a receptor and display the result in Nanome.

### Installation

```sh
$ pip install nanome-docking
```

### Usage

To start the plugin:

```sh
$ nanome-docking smina -a plugin_server_address
```

Runs only on Linux

OR

```sh
$ nanome-docking autodock4 -a plugin_server_address
```

Runs on Windows, and needs Autodock4 to be installed on the computer and in the PATH variable.

In Nanome:

- Activate Plugin, Docking window will automatically open
- Select a receptor
- Click on "Ligand", and select ligands to dock
- If using Smina, click on "Site", and select which molecule should be used to define the docking site
- Choose exhaustiveness, number of results and size of the box to generate around the site molecule (Smina only)
- If visual scoring is turned on, atom size and labels will indicate each atom's contribution to the ligand's score
- Click on Run

### License

MIT