# SOFA Tearing Plugin
This plugin is for SOFA to implement the tearing for Triangular Meshes. A tear threshold value is considered in plugin implementation to obtain how much force needs to cause the mesh to tear. The position and orientation of the tear calculate in each step and then mesh topology change and tear take place.<br /><br />
This plugin is compatible with **SOFA v21.06**.<br />
You can find documentation and tutorials of how to use SOFA in https://www.sofa-framework.org/.<br />


## How to compile it?
There is a video that explains how to compile an external plugin [here](https://youtu.be/46E215871e8).<br /><br />
Below there are steps to build this plugin.
#### 1. Compile SOFA
Download and compile sofa source and make sure **runSofa** is working.<br />
You can use this [toturial](https://www.sofa-framework.org/community/doc/getting-started/build/windows/) to compile SOFA source.
#### 2. Making directory for external plugins
You want compile and use one or more external plugins it is preferred to create one specific repository outside SOFA where you can checkout all these external plugins. This structure is preferred since it will allow a clean organization of external plugins in one single repository. Let’s note the path to this repository */ext_plugin_repo/*.<br /><br />
In this directory, the structure is:
- */ext_plugin_repo/*
  - */SofaTearing*
#### 3. Clone plugin source code
Clone this repository into SofaTearing directory that you created in previous step.
#### 4. CMakeList of the repository
In order to handle this repository as one single set of external plugins, you need to write a short CMakeList.txt file in */ext_plugin_repo/* as follows:<br />
```
cmake_minimum_required(VERSION 3.12)

find_package(SofaFramework)

sofa_add_plugin(SofaTearing/  SofaTearing)
```
#### 5. CMake option in SOFA
To compile all the external plugins located in this repository, all you need to do is to set the path to this repository */ext_plugin_repo/* in the CMake variable: SOFA_EXTERNAL_DIRECTORIES.

This will directly configure and allow to compile all specified plugins from SOFA.

## How to use it
There is an example of using this plugin in your scene file in the example directory.