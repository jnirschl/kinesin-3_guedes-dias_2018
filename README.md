# Supplemental code for "Kinesin-3 responds to local microtubule dynamics to target synaptic cargo delivery to the presynapse"

[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](TODO--add_url/blob/master/LICENSE) [![DOI](https://zenodo.org/badge/doi/00.0000/zenodo.00000.svg)](http://dx.doi.org/00.0000/guedesdias.00000)

------------------
## Contents
* [Installation](#installation)
* [Files](#files)
	* [ensembleRecruitment](ensembleRecruitment): A script to segment microtubules and quantify fluorescence intensity of objects on microtubules vs. background to measure microtubule binding in TIRF microscopy.
	* [kymoSuite](kymoSuite): A set of MALTAB (R) functions with GUI for manual kymograph tracing of motility events and semi-automated analysis.
	* [nearestPresynapticDistance](nearestPresynapticDistance): A script to measure the observed nearest neighbor distance from one set of *n* points (e.g., EB3 comet initiation/ termination) to a set of reference *p* points (e.g., presynapses) along a 1-dimensional line *l*. For each set of observations, it performs a Monte Carlo simulation of *n* points distributed with a uniform random probability along *l* and computes the nearest neighbor distance to a fixed reference set *p*.
* [License and citation](#license-and-citation)
* [Acknowledgements](#acknowledgements)

------------------
## Installation
These scripts use the MATLAB (R) programming language and environment. To begin using these scripts, simply add the files to the MATLAB path using the "Set Path" dialog box or the `addpath` function. Previous experience with MATLAB (R) or another programming language is required for efficient use of these scripts.

### System requirements

| System requirements |                   |
| ----------          | :----------:      |
| Operating system    | > Windows 8.1     |
|                     | > Ubuntu 16.04    |
| Software            | MATLAB (R)        |
| Processor           | 2 GHz             |
| Memory (RAM)        | 2 GB              |
| Hard-drive space    | 50 MB             |

**Minimal system requirements for using these scripts.** This does not include hard-drive space required for additional software packages or dependencies.

### Dependencies
The following MATLAB (R) Toolboxes are required: Image Processing Toolbox, Statistics and Machine Learning Toolbox, and Signal Processing Toolbox. Additional dependencies include [uTrack (Danuser lab)](https://www.utsouthwestern.edu/labs/danuser/software/#utrack_anc), [Bioformats](https://docs.openmicroscopy.org/bio-formats/5.9.2/users/matlab/), and MATLAB Central files where listed.


------------------
## License and Citation
These scripts are released under the [MIT License](https://opensource.org/licenses/MIT).

If you use use this code in your research, please cite:
Guedes-Dias P, Nirschl JJ, Abreu N, Tokito MK, Magiera MM, Janke C, and Holzbaur ELF. Kinesin-3 responds to local microtubule dynamics to target synaptic cargo delivery to the presynapse. Current Biology. 2018.

Bibtex formatted reference:
```text
@article{GuedesDias2018,
    Author={Guedes-Dias P, Nirschl JJ, Abreu N, Tokito MK, Magiera MM, Janke C, and Holzbaur ELF.},
    Journal={Current Biology},
    Title={Kinesin-3 responds to local microtubule dynamics to target synaptic cargo delivery to the presynapse},
    Year{2018},
}
```

------------------
## Acknowledgements
This work was funded by NIH Grants R01GM048661 (to E.L.F.H.) and F30NS092227 (to J.J.N.). This study has received support under the program “Investissements d’Avenir” launched by the French Government and implemented by ANR with the references ANR-10-LBX-0038, ANR-10-IDEX-0001-02 PSL. The work of M.M. is supported by the Fondation pour la Recherche Medicale (FRM) grant DEQ20170336756 and the the Fondation Vaincre Alzheimer grant FR-16055p. P.G-D. was supported by the NSF Science and Technology Center for Engineering MechanoBiology (Grant CMMI-1548571).