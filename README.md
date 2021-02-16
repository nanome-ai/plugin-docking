# Nanome - Docking

This Nanome plugin interfaces with a variety of docking softwares to dock ligands to a receptor and display the result in Nanome.

### Installation



### Docker Usage

To run in a Docker container:

```sh
$ cd docker
$ ./build.sh

Note that if you get the error "no permission to read from '.../plugin-docking/nanome_docking/smina' ", 
change the permission using:
$ chmod a+r ../nanome_docking/smina
and try the previous step( $./build.sh ) again:

$ ./deploy.sh -a <plugin_server_address> smina [args]
```




To view (and follow the logs of) a container:
```sh
$ docker logs --follow docking
```

To stop the container:
```sh
$ docker stop docking
```

To view the running containers:
```sh
$ docker ps
```


### Install from Pip or Pip3


```sh
$ pip install nanome-docking
```


To start the plugin:

```sh
$ nanome-docking smina -a <plugin_server_address>
```

Runs only on Linux

OR

```sh
$ nanome-docking autodock4 -a <plugin_server_address>
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


### Run in development

```
python3 run.py -r -a plugins.nanome.ai -p 9999
```



### License

MIT
