# Notebook for Network and Topological Analysis in Neuroscience
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![GLP-3.0 License](https://img.shields.io/github/license/multinetlab-amsterdam/network_TDA_tutorial)](https://img.shields.io/github/license/multinetlab-amsterdam/network_TDA_tutorial)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)


[![Stars](https://img.shields.io/github/stars/multinetlab-amsterdam/network_TDA_tutorial?style=social)](https://img.shields.io/github/stars/multinetlab-amsterdam/network_TDA_tutorial?style=social)
[![Watchers](https://img.shields.io/github/watchers/multinetlab-amsterdam/network_TDA_tutorial?style=social)](https://img.shields.io/github/watchers/multinetlab-amsterdam/network_TDA_tutorial?style=social)

[![DOI](https://zenodo.org/badge/321418758.svg)](https://zenodo.org/badge/latestdoi/321418758)



Authors: Eduarda Centeno  & Fernando Santos 

Contact information: <e.centeno@amsterdamumc.nl> or <f.nobregasantos@amsterdamumc.nl>

-------------------------------------------------------------------------

## Table of contents
1. [General information](#general_information)
2. [Requirements](#requirements)
3. [How to install](#how_to_install)
4. [Notes](#notes)
5. [Acknowledgements](#acknowledgements)

### <a id='general_information'></a> General information:

 The primary purpose of this project is to facilitate the computation of Network Neuroscience metrics using Graph Theory and Topological Data Analysis.
 
<br>

> This repository is a supplement to the preprint [A hands-on tutorial on network and topological neuroscience](https://doi.org/10.1101/2021.02.15.431255) and the [final version published at Brain Structure and Function](https://doi.org/10.1007/s00429-021-02435-0).

<br>

#### We divided this notebook into two parts:

1. The first part contains the standard computations of network and TDA metrics and visualizations. 
2. The second part is dedicated to the 3D visualizations developed by our group.


```
                  A video created with our 3D brain plots!
```
![](./Figures/filtration.gif)


------------------------------------------------------------------------

### <a id='requirements'></a> Requirements:
Here we will describe the core packages that will be necessary, but due to dependency compatibilities, we have provided a requirements.txt with all packages needed to be installed in a new Anaconda environment. 

    - Python: 3.x
    - Numpy: 1.18.5
    - Matplotlib: 3.3.2
    - Meshio: 4.0.16 --- https://pypi.org/project/meshio/
    - Seaborn: 0.11.0
    - Pandas: 1.1.3
    - Networkx: 2.4
    - Nxviz: 0.6.2  --- For CircosPlot to work fully, we recommend installing through https://github.com/eduardacenteno/nxviz/
    - Community (python-louvain): 0.13 --- https://python-louvain.readthedocs.io/en/latest/api.html
    - Gudhi: 3.4.1 --- http://gudhi.gforge.inria.fr/
    - Plotly: 4.6.0
    - Scikit-learn: 0.24.1
    - Scipy: 1.5.0
    - JupyterLab: 1.2.0 --- This is very important; otherwise, plotly will not work as we intended. (https://plotly.com/python/getting-started/)

-------------------------------------------------------------------------

### <a id='how_to_install'></a>How to install:
We recommend creating a new environment in Anaconda dedicated for the use of these notebooks (Primarly tested on Windows/Linux).

1. Create a new Anaconda environment with the correct python version (in Anaconda prompt or navigator)

2. Activate the new environment in the command line (Anaconda prompt)

```
conda activate envname
```

3. Change to the notebook's directory

```
cd path\to\notebookfolder
```

4. Install packages using pip in Anaconda prompt

** with environment-specific python.exe (Windows)
```
path\to\anaconda3\envs\envname\python.exe -m pip install -r requirements.txt
```
** MacOS users do not need to give the path to the environment's python.exe
```
pip install -r requirements.txt
```

5. Add jupyter-plotly labextension (key for 3D visualization)

```
jupyter labextension install jupyterlab-plotly 
```
<br>

__Troubleshooting__:

1. Permission error and suggestion to use --user (possibly because the user did not use environment-specific python.exe)
    
    > Try opening the prompt command as admin.
    
    > Using *pip install anaconda* before installing packages will potentially solve the issue.
    
2. Jupyter Lab is asking for Node.js 5+. 

    > Using *conda install nodejs* will potentially solve the issue.
    
<br>

#### Web-based options:
| Nbviewer | Jupyter Notebook | Jupyter Lab\*| HTML |
| ---      | --               | ---         | ---  |
| [1-network_analysis.ipynb](https://nbviewer.jupyter.org/github/multinetlab-amsterdam/network_TDA_tutorial/blob/main/1-network_analysis.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/multinetlab-amsterdam/network_TDA_tutorial/HEAD?filepath=1-network_analysis.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/multinetlab-amsterdam/network_TDA_tutorial/HEAD?urlpath=lab/tree/1-network_analysis.ipynb) | [HTML](https://htmlpreview.github.io/?https://github.com/multinetlab-amsterdam/network_TDA_tutorial/blob/main/html/1-network_analysis.html) |
| [2-visualization_3d.ipynb](https://nbviewer.jupyter.org/github/multinetlab-amsterdam/network_TDA_tutorial/blob/main/2-visualization_3d.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/multinetlab-amsterdam/network_TDA_tutorial/HEAD?filepath=2-visualization_3d.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/multinetlab-amsterdam/network_TDA_tutorial/HEAD?urlpath=lab/tree/2-visualization_3d.ipynb) | [HTML](https://htmlpreview.github.io/?https://github.com/multinetlab-amsterdam/network_TDA_tutorial/blob/main/html/2-visualization_3d.html) |

#### \* Binder is not building correctly plotly extensions in Jupyter Lab; therefore, 3D visualizations in the second notebook can only be seen using the jupyter notebook option!
-------------------------------------------------------------------------

### <a id='Notes'></a>Notes:

The jupyter notebooks can throw some warnings due to package updates and resulting deprecations. It is possible to use the following code lines to ignore these warnings:

```
import warnings

warnings.filterwarnings('ignore') 
```

-------------------------------------------------------------------------

### <a id='acknowledgements'></a>Acknowledgements:
The 1000_Functional_Connectomes dataset was downloaded from the [The UCLA multimodal connectivity database](http://umcd.humanconnectomeproject.org/). 

Brown JA, Rudie JD, Bandrowski A, Van Horn JD, Bookheimer SY. The UCLA multimodal connectivity database: a web-based platform for brain connectivity matrix sharing and analysis. Frontiers in neuroinformatics. 2012 Nov 28;6:28.  (http://dx.doi.org/10.3389/fninf.2012.00028) 

Biswal BB, Mennes M, Zuo XN, Gohel S, Kelly C, Smith SM, Beckmann CF, Adelstein JS, Buckner RL, Colcombe S, Dogonowski AM. Toward discovery science of human brain function. Proceedings of the National Academy of Sciences. 2010 Mar 9;107(10):4734-9. [Freely available dataset](http://fcon_1000.projects.nitrc.org/)

Rule A, Birmingham A, Zuniga C, Altintas I, Huang SC, Knight R, Moshiri N, Nguyen MH, Rosenthal SB, Pérez F, Rose PW. Ten simple rules for writing and sharing computational analyses in Jupyter Notebooks. https://github.com/jupyter-guide/ten-rules-jupyter
