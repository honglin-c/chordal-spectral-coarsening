# Chordal Decomposition for Spectral Coarsening

This is the MATLAB implementation of "Chordal Decomposition for Spectral Coarsening" [Chen et al. 2020]. The only dependency is the gptoolbox: `https://github.com/alecjacobson/gptoolbox`.

* main_specCoarse.m: demo for comparison with [Liu et al. 2019]  on spectral coarsening of cotangent Laplacian

* main_qslim_cotan.m: demo for optimizing cotangent laplacian starting from the coarsened Laplacian of [Garland & Heckbert 1997]

* main_vol2surf.m: demo for volume-to-surface approximation of the Laplacian

* main_qslim_anisotropic.m: demo for optimizing anisotropic Laplacian starting from the coarsened Laplacian of [Garland & Heckbert 1997]

This program has been tested on MATLAB R2019b and R2020a.

## Bibtex

```
@article{Chen:ChordalSpecCoarsen:2020,
  title = {Chordal Decomposition for Spectral Coarsening},
  author = {Honglin Chen and Hsueh-Ti Derek Liu and Alec Jacobson and David I.W. Levin},
  year = {2020},
  issue_date = {December 2020}, 
  publisher = {Association for Computing Machinery}, 
  volume = {39}, 
  number = {6}, 
  issn = {0730-0301},
  journal = {ACM Trans. Graph.}, 
}
```

