# Automatic-Saltykov
## Scope
When dealing with polycrystals, such as metals or minerals, the microstructures are usually described from 2D sections 
(e.g. thanks to optical microscopy). Thus, the grain size distribution is chatacterized by the apparent grains sizes only. 
Unfolding the apparent grain size distribution in order to evaluate the actual 3D distribution is usually referred as the corpuscle problem.

Assuming that the grains are equiaxed enough to treat them like spheres, the Scheil-Schwartz-Saltykov method (called Saltykov method for 
short) can be used [[1]](#1). The Saltykov method is an iterative method, working on a finite histogram and based on the Wicksell's 
equation [[2]](#2). Nevertheless, this method comes with a couple of options (e.g. the number of classes to be used), and there is no 
general rule about to choose those options. 

This programs aims to automatically find the best options for the Saltykov method. Here, the results are considered as "the best" if the
*refolded* distribution minimizes the Cramer-von Mises good-of-fit test [[3]](#3).

## Usage
Given a sample of apparent radii ``R``, the following command performs the Saltykov method with automatic parameters and plots the 
unfolded histogram:

    autoSaltykov(R)
    
Look at the documentation for further details about the input and output options:

    help autoSaltykov(R)
    
## Cite this work
When using the present functions in your research papers, please cite reference [[3]](#3). Alternatively, use the following BibTeX 
entry:

    @article{ImageAnalStereol2133,
	    author = {Dorian Depriester and R\'{e}gis Kubler},
	    title = {Resolution of the Wicksell's equation by Minimum Distance Estimation},
	    journal = {Image Analysis & Stereology},
	    volume = {38},
	    number = {3},
	    year = {2019},
	    issn = {1854-5165},	
      pages = {213--226},	
      doi = {10.5566/ias.2133},
    }

## References
<a id="1">[1]</a> Higgins, M.D (2000). Measurement of crystal size distributions. American Mineralogist, 85 (9): 1105–1116. doi:[10.2138/am-2000-8-901](https://doi.org/10.2138/am-2000-8-901)

<a id="2">[2]</a> Wicksell, S. (1925). The Corpuscle Problem: A Mathematical Study of a Biometric Problem. Biometrika, 17(1/2), 84-99. doi:[10.2307/2332027](https://doi.org/10.2307/2332027)

<a id="3">[3]</a> Depriester, D., & Kubler, R. (2019). Resolution of the Wicksell's equation by Minimum Distance Estimation. Image Analysis & Stereology, 38(3), 213-226. doi:[10.5566/ias.2133](https://doi.org/10.5566/ias.2133)

