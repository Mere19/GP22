# Assignment 4

Name: 'Di Zhuang'

Legi-Nr: '21958772'

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Mandatory Tasks

1) Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off, Octo_cut2.obj)

2) Several examples of the distortion visualizations.


## Reports
### (mandatory-1) parameterization and checkerboard texture models
#### cathead
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/cathead_uni_texture.png" width="300">| <img align="center"  src="./res/cathead_uni_param.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/cathead_cot_texture.png" width="300">| <img align="center"  src="./res/cathead_cot_param.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/cathead_lscm_texture.png" width="300">| <img align="center"  src="./res/cathead_lscm_param.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/cathead_arap_texture.png" width="300">| <img align="center"  src="./res/cathead_arap_param.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/cathead_lscm_texture_free.png" width="300">| <img align="center"  src="./res/cathead_lscm_param_free.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/cathead_arap_texture_free.png" width="300">| <img align="center"  src="./res/cathead_arap_param_free.png" width="300"> |

#### hemisphere
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/hemi_uni_texture.png" width="300">| <img align="center"  src="./res/hemi_uni_param.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/hemi_cot_texture.png" width="300">| <img align="center"  src="./res/hemi_cot_param.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/hemi_lscm_texture.png" width="300">| <img align="center"  src="./res/hemi_lscm_param.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/hemi_arap_texture.png" width="300">| <img align="center"  src="./res/hemi_arap_param.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/hemi_lscm_texture_free.png" width="300">| <img align="center"  src="./res/hemi_lscm_param_free.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/hemi_arap_texture_free.png" width="300">| <img align="center"  src="./res/hemi_arap_param_free.png" width="300"> |

#### hemisphere_non_convex_boundary
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/heminonconv_uni_texture.png" width="300">| <img align="center"  src="./res/heminonconv_uni_param.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/heminonconv_cot_texture.png" width="300">| <img align="center"  src="./res/heminonconv_cot_param.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/heminonconv_lscm_texture.png" width="300">| <img align="center"  src="./res/heminonconv_lscm_param.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/heminonconv_arap_texture.png" width="300">| <img align="center"  src="./res/heminonconv_arap_param.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/heminonconv_lscm_texture_free.png" width="300">| <img align="center"  src="./res/heminonconv_lscm_param_free.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/heminonconv_arap_texture_free.png" width="300">| <img align="center"  src="./res/heminonconv_arap_param_free.png" width="300"> |

#### Octo_cut2
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/octo_uni_texture.png" width="300">| <img align="center"  src="./res/octo_uni_param.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/octo_cot_texture.png" width="300">| <img align="center"  src="./res/octo_cot_param.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/octo_lscm_texture.png" width="300">| <img align="center"  src="./res/octo_lscm_param.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/octo_arap_texture.png" width="300">| <img align="center"  src="./res/octo_arap_param.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/octo_lscm_texture_free.png" width="300">| <img align="center"  src="./res/octo_lscm_param_free.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/octo_arap_texture_free.png" width="300">| <img align="center"  src="./res/octo_arap_param_free.png" width="300"> |


### (mandatory-2) distortion visualization
#### cathead
| mtd \ metric      | Conformal (angle) |    Authalic (area)  |  Isometric  (length)    |
| :--------------:  | ----------------- | ------------------- | ----------------------- |
| LSCM (free)       |<img align="center" src="./res/cathead_lscm_angle.png" width="300">| <img align="center"  src="./res/cathead_lscm_area.png" width="300"> | <img align="center"  src="./res/cathead_lscm_length.png" width="300"> |
| ARAP (free) |<img align="center" src="./res/cathead_arap_angle.png" width="300">| <img align="center"  src="./res/cathead_arap_area.png" width="300"> |<img align="center"  src="./res/cathead_arap_length.png" width="300"> |


#### hemisphere
(copy the above table to format)
| mtd \ metric      | Conformal (angle) |    Authalic (area)  |  Isometric  (length)    |
| :--------------:  | ----------------- | ------------------- | ----------------------- |
| LSCM (free)       |<img align="center" src="./res/hemi_lscm_angle.png" width="300">| <img align="center"  src="./res/hemi_lscm_area.png" width="300"> | <img align="center"  src="./res/hemi_lscm_length.png" width="300"> |
| ARAP (free) |<img align="center" src="./res/hemi_arap_angle.png" width="300">| <img align="center"  src="./res/hemi_arap_area.png" width="300"> |<img align="center"  src="./res/hemi_arap_length.png" width="300"> |
