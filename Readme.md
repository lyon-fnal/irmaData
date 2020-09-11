# irmaData

Adam Lyon (September 2020)

Get data for IRMA analysis into HDF5 and analyze at an HPC facility. 

- [irmaData](#irmadata)
	- [1 Introduction](#1-introduction)
	- [2 Development](#2-development)
		- [2.1 Local host with no dependencies](#21-local-host-with-no-dependencies)
			- [2.2.1 Code tours](#221-code-tours)
		- [2.2 Local host with dependencies](#22-local-host-with-dependencies)
		- [2.3 Remote SL7 node](#23-remote-sl7-node)

## 1 Introduction

This repository contains code for converting Muon g-2 data into the HDF5 format and understanding what data is required. Development is in the `art-studio` environment within the `project` directory. See the `setup.sh` script within. 

In `project/src` there are three modules...

- `ratioeast_module.cc` is a fairly recent version of Sudeshna's module for making IRMA histograms.
- `james_module.cc` is James' version of the histogram making module.
- `irmaData.cc` is the module that creates an HDF5 file with the data necessary to make histograms (determined by looking at James' module).

If you want to make HDF5 files, you'll run `irmaData.cc` with the `irmaData.fcl` FHICL file. 

## 2 Development

You can develop within this repository several ways (and see sections below for more details),

- Local host with no dependencies: Use Microsoft Visual Studio Code (VSCode) and develop in a sandbox on your host (you need Docker and VSCode installed - that's it)
- Local host with dependencies: Clone to your local host and, with VSCode, attach to a `devenv` container you've already installed and have up
- Remote SL7 node: Clone to a Scientific Linux 7 machine with CVMFS (e.g. `gm2testgpvm01`)

### 2.1 Local host with no dependencies

Use this choice if you want to try out this repository without installing it permanently and you don't have an installation of [devenv](https://github.com/lyon-fnal/devenv). The only things you need on your machine are Microsoft Visual Code (VSCode) and Docker.

- Install Docker for [Mac](https://hub.docker.com/editions/community/docker-ce-desktop-mac/) or [Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows/). Get the stable version. 
- Install [Microsoft Visual Studio Code](https://code.visualstudio.com) (VSCode) (choose the stable version)
- Start up VSCode
- Click on the gear icon in the lower left corner and select Extensions. Search for and install the `Remote - Containers` extension (you may want the `Remote-SSH` one too).
- From the Command Pallette (Command-Shift-P on Mac, Control-Shift-P on Windows), choose "Remote-Containers: Open Repository in Container". The repository name is `lyon-fnal/irmaData`
- Choose `Create a new volume`
- It may take awhile to start if docker is downloading the image
- When things start up, mount the CVMFS volumes...
  - Bring up a terminal with Control-` (back quote)
  - Type `/usr/local/bin/start_cvmfs.sh`
- Load the environment by clicking on `gm2-sl7-v9-irmaData` on the bottom status bar and selecting that same "kit" from the list
- Configure CMake by bringing up the Command Palette and choosing `CMake: Configure`. It may take awhile to run because you are populating your CVMFS cache
- Build the code by bringing up the Command Palette and choosing `CMake: Build`. Building may take awhile
- In order to run, activate a terminal with Control-` (back quote) and do the following...
  - `source setup.sh`
  - `mkdir try ; cd try`
  - Get a kerberos ticket with `kinit`
  - Make mount point for `/pnfs` with `mkdir /pnfs`
  - Mount `/pnfs` with `sshfs`. For example, (replace `you` and `machine` accordingly)...

```bash
sshfs -o reconnect,ServerAliveInterval=5,ServerAliveCountMax=2 you@machine:/pnfs /pnfs
```

Run with (for example)

```bash
gm2 -c irmaData.fcl -n 10 /pnfs/GM2/daq/run2/offline/gm2_5111C/runs_25000/25024/gm2offline_final_31680117_25024.00365.root
```

It will take awhile to run for the first time as you still have to populate your CVMFS cache.
You can look at the structure of the resulting `irma.h5` file with

```bash
h5ls-shared irma.h5/ReconEastClusters
h5ls-shared -v irma.h5/ReconEastClusters/event
```

#### 2.2.1 Code tours

The *Code Tour* VSCode extension is automatically installed when you launch remote development into the container. Code tours are curated  and annotated "walks" through code. It's a nice way to learn how code works. There are two code tours in this repository. To view them, bring up the file side-bar (should be the default, but if not, click the top document looking icon on the left of the window). You should see a heading towards the bottom left titled `CODETOUR`. Click on the little arrow next to the heading. The two code tours are,

- James Histogram Maker - A tour of James' analyzer module to learn how it writes out the histograms.
- Make HDF5 File - A tour of the analyzer module that makes HDF5 files.

### 2.2 Local host with dependencies

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


### 2.3 Remote SL7 node

Log into the SL7 node. It needs to have CVMFS installed. `cd` to your development area. 

Simply clone the repository the normal way.  Before you source the `setup.sh` script, you will need to set the `WORKSPACE_FOLDER` environment variable to the path to the top directory in the repository. 
Then you can set up the environment with source `setup.sh`, build with `studio build` and run.

## 3 Status

I ran over era 2E (20K files).