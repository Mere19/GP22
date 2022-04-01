# Assignment 3

Name: 'Your real name'

Legi-Nr: 'Your legi number'

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Mandatory Tasks
1) Show the screenshots of mesh **torus.obj** and the mesh **Julius.obj** from the data folder shaded with the 5 different normals. For the PCA and quadratic fitted normals show them each with both k=1 and k=2. Note that here k means k-ring neighbors.

2) Report your matching between curvature expressions and the figures. Provide a short motivation for your matchings.

3) Show screenshots of mesh **bumpy-cube.obj** from the data folder coloured according to the 4 discrete curvature measures.

4) Report your findings on the differences between implicit and explicit Laplacian smoothing, the differences between uniform vs cotangent weights and report the parameters you used and how changing them affects the result. Also show screenshots of both the implicit and explicit Laplacian smoothing results for mesh **bunny_noise.obj** from the data folder.

5) Report your findings on the comparison with results you obtained with the Laplacian smoothing vs with bilateral smoothing. Show screenshots of bilateral smoothing applied to mesh **bunny_noise.obj** from the data folder.


## Reports
### 1 - Shading w.r.t. different normals

**Use the mesh torus.obj and the mesh Julius.obj**

**Use the built-in function igl::per_vertex_normals() to orient your normals consistently**

| normals        | torus.obj                  | Julius.obj                 |
| :-----------:  | ------------------------------------- |------------------------------------- |
| standard       |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| area-weighted  |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| mean-curvature |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| PCA (k=1)      |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| PCA (k=2)      |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| quadratic (k=1)|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| quadratic (k=2) |<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |

### 2 - Matching curvature expressions and the figures
| expression   |  Your answer |  Your explanation   |
|--------------|--------------|------------------|
| k1           | a/b/c/d      | xxxxx            |
| k2           | a/b/c/d      | xxxxx            |
| k3           | a/b/c/d      | xxxxx            |
| k4           | a/b/c/d      | xxxxx            |


### 3 - Visualize curvatures

**Use the mesh bumpy-cube.obj**

| Min Curvature                         |  Max Curvature                       |
| ------------------------------------- |------------------------------------- |
|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |
| Mean Curvature                        |  Gaussian Curvature                  |
|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |


### 4 - Implicit v.s. explicit Laplacian Smoothing

**Use the mesh bunny_noise.obj**

**Try different laplacian matrices, step sizes and iterations**

| Input  |  Implicit (your params01)    |  Implicit (your params02)          | Implicit (your params03)          |
| -------|----------------------------- |------------------------------------|---------------------------------- |
|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |<img align="center"  src="./res/placeholder.png" width="300"> |<img align="center"  src="./res/placeholder.png" width="300"> |

Your observations:


| Input  |  Explicit (your params01)    |  Explicit (your params02)          | Explicit (your params03)          |
| -------|----------------------------- |------------------------------------|---------------------------------- |
|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |<img align="center"  src="./res/placeholder.png" width="300"> |<img align="center"  src="./res/placeholder.png" width="300"> |

Your observations:

### 5 - Laplacian v.s. bilateral smoothing

**Use the mesh bunny_noise.obj**

| Input                                 |  Laplacian Smoothing                 |  Bilateral Smoothing                 |
| ------------------------------------- |------------------------------------- |------------------------------------- |
|<img align="center" src="./res/placeholder.png" width="300">| <img align="center"  src="./res/placeholder.png" width="300"> |<img align="center"  src="./res/placeholder.png" width="300"> |

Your observations:
