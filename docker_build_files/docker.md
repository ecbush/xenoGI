# xenoGI: Code for reconstructing genome evolution in clades of microbes

Source code is available at https://github.com/ecbush/xenoGI

This is a basic introduction to how to run xenoGI from within docker. It's intended for those with little familiarty with docker.

## Downloading docker

The first step is to install docker on your machine. That's available here: https://www.docker.com/products/docker-desktop/

## Mounting the docker image
### Preparing docker settings

Once you've started the docker desktop application, you can go to the gear (upper right). Then under "Resources" you can specify how many of your machine's CPUs and how much RAM to allocate to docker.

In general, more CPUs are better, because xenoGI is able to use multiple threads. Once you have set a number of processors, you can then enter that same number in the the xenoGI ``params.py`` file. If you have 8 CPUs allocated, then you would seet ``numProcesses`` to 8 in that file.

RAM usage depends on how many strains you have. For the example data set, 6 GB should be sufficient. If you later use a larger data set and don't allocate enough memory, then you'll get an error message saying the process has been killed.

### To mount the image as a container, perform the following steps:
Pull the docker image from docker hub.

If you are running on an amd64 machine (most people) then do this:

          $ docker pull ecbush/xenogi:3.1.0_amd64

If you are running on an arm64 machine (ie an M1 mac) then do this:

          $ docker pull ecbush/xenogi:3.1.0_arm64

You can then create a directory where you're indending to have xenoGI run. Let's say you have put a copy of the ``example/`` directory from the github repository in /home/data/example. You can run a container (instance) of the image you just got, and mount that example directory into its file system. Like this:

          $ docker run -it --name myContainer -v /home/data/example:/example ecbush/xenogi:3.1.0_amd64 bash

Note that if you were on an M1 mac, then you'd substitute arm64 for amd64. Also note that you must provide the __ABSOLUTE PATH__ to your working directory.

## Running xenoGI

If you have successfully mounted the image, then you should see something like this:

          root@fe46a3c61f0b:/#

You should then be able to cd into the example directory and run xenoGI:

          root@fe46a3c61f0b:/# cd example
          root@fe46a3c61f0b:/example# xenoGI params.py runAll

After the run completes, you can do interactiveAnalysis etc as described in the README and TUTORIAL on the github repository. When you are done, you can quit out of the docker container by typing ``exit()`` or ``<ctrl> d``. And you can remove the container via the docker desktop app. (or by running ``docker container rm myContainer`` at the command line) 

The xenoGI directory you have worked on will remain, and can be accessed on your computer's file system. (e.g. in our example this would be at ``/home/data/example``).