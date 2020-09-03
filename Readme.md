# irmaData
Get data for IRMA analysis

## Development
You can develop this repository several ways,

- Local host with no dependencies: Use Microsoft Visual Studio Code (VSCode) and develop in a sandbox on your host (you need Docker and VSCode installed - that's it)
- Local host with dependencies: Clone to your local host and, with VSCode, attach to a `devenv` container you've already installed and have up
- Remote SL7 node: Clone to a Scientific Linux 7 machine with CVMFS (e.g. `gm2testgpvm01`)

### Local host with no dependencies

Use this choice if you want to try out this repository without installing it permanently and you don't have an installation of [devenv](https://github.com/lyon-fnal/devenv). The only things you need on your machine are Microsoft Visual Code (VSCode) and Docker.

- Start up VSCode
- From the Command Pallette, choose "Remote-Containers: Open Repository in Container". The repository name is `lyon-fnal/irmaData`
- Choose `Create a new volume`

It may take awhile to start if docker is downloading the image. 

Once in the container, you should wait a minute or two for CVMFS to mount. You can check with `ls /cvmfs/gm2.opensciencegrid.org/`. Once that is done, you should select the `gm2-sl7-v9-irmaData` kit and then Configure. If you get an error when Configuring, you may need to choose `CMake: Delete Cache and Reconfigure` from the command pallette. Configuring again usually works. Then build. 

To run, you'll need to open the integrated terminal and source the `setup.sh` script first. 

Note that things will be slow the first time you configure, build, and run due to populating the CVMFS cache. 

### Local host with dependencies

Start up the `devenv` container (see [devenv](https://github.com/lyon-fnal/devenv)). If you aren't going to use VSCode, then simply exec into it and clone this repository. Follow the instructions for *Remote SL7 node*. 

If you want to use VSCode, then do the following...

- Clone this repository to your local disk
- Start VSCode
- Bring up the command pallette and choose "Remote-Containers: Attach to running container" and choose the container
- Bring up the command pallette and do "Remote-Containers: Open Named Container Configuration File"
- Add the following within the outermost braces...

```json
"workspaceFolder": "/Users/lyon/Development/gm2/irmaData",
	"remoteEnv": {
		"WORKSPACE_FOLDER":"/Users/lyon/Development/gm2/irmaData",
		"USER":"lyon"
	},
	"extensions": [
		"ms-vscode.cmake-tools",
		"ms-vscode.cpptools",
		"twxs.cmake"
    ]
```

where you replace the paths to point to where you cloned this repository and change `USER` accordingly. 

You will need to close your connection to the container and reattach for this to go into effect. 

Once you've done that, you can use VSCode. Select the `gm2-sl7-v9-irmaData` kit, Configure, and Build. 


### Remote SL7 node

Log into the SL7 node. It needs to have CVMFS installed. `cd` to your development area. 

Simply clone the repository the normal way.  Before you source the `setup.sh` script, you will need to set the `WORKSPACE_FOLDER` environment variable to the path to the top directory in the repository. 
Then you can set up the environment with source `setup.sh`, build with `studio build` and run.
