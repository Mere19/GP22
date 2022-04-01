# Assignment 2
Name: 'Di Zhuang'

Legi-Nr: '21958772'

**The ouptut .off files are in assignment2/res folder**

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the .off format for every point-cloud dataset provided (assignment2/res).

4) Theory question (1): Save your notes to assignment2/res and link to them from within the report.

5) Theory question (2): Save your notes to assignment2/res and link to them from within the report.

6) Show screenshots comparing the 'hound.off' computed with normal-based reconstruction to the point-based reconstruction of the previous task.

7) Compare your MLS reconstruction results to the surfaces obtained with RIMLS and Screened Poisson Reconstruction, and inspect the differences. Report your findings.

8) Show screenshorts of your method for estimating normals based on Principal Component Analysis. Compare the results with the reconstruction obtained using the provided normals.

## Reports
### 1 - Visualization of the 'cat.off' point cloud with constrained points
| model name  | view 01             | view 02            |
| :---------: | ------------------- | ------------------ |
| cat         |<img align="center" src="./res/Q1_v1.png" width="300">| <img align="center"  src="./res/Q1_v2.png" width="300"> |

The difference between spatial-index method and no-spatial-index method is not obvious for cat.off. Instead, I tested these two methods on a denser point cloud horse.off, with resolution of the grids equal to epsilon. This gives 0.16 seconds for spatial-index method and 4 seconds for non-spatial-index method.

### 2 - Grid with nodes colored w.r.t. the implicit function (using the non-axis-aligned grid as described in Section 2.3) 
| model name  | view 01             | view 02            |
| :---------: | ------------------- | ------------------ |
| cat (polyDegree = 0, resolution = 30, wendlandRadius = 100)         |<img align="center" src="./res/Q2_cat_v1.png" width="300">| <img align="center"  src="./res/Q2_cat_v2.png" width="300"> |
| luigi (polyDegree = 0, resolution = 40, wendlandRadius = 10)      |<img align="center" src="./res/Q2_luigi_v1.png" width="300">| <img align="center"  src="./res/Q2_luigi_v2.png" width="300"> |

### 3 - Reconstructed surfaces
**Please also save the reconstructed shape (.off) in the results folder**

| sphere                                                     | cat                                                           |
| :--------------------------------------------------------: |  :----------------------------------------------------------: | 
|<img align="center" src="./res/Q3_sphere.png" width="300">| <img align="center"  src="./res/Q3_cat.png" width="300"> |
| luigi                                                      | horse                                                         |
|<img align="center" src="./res/Q3_luigi.png" width="300">| <img align="center"  src="./res/Q3_hors.png" width="300"> |
| hound                                                     |                                                               |
|<img align="center" src="./res/Q3_hound.png" width="300">|                                                               |


**Please summarize your observations of what happens when you change the following parameters. Please feel free to add screenshots (by creating a separate table) to illustrate your explanation.**

| params                   | Your Observation    | 
| :---------------------:  | ------------------- |
| grid resolution          |   As grid resolution increases, the reconstructed mesh becomes more smooth.   |
| Wendland function radius |   As the Wendland function radius increases, the continuity of the reconstructed mesh becomes better.   |
| polynomial degree        |  As the polynomial degree increases, the details of the surface are better approximated. However, increasing the polynomial degree may also lead to overfitting, resulting in noise in the reconstructed mesh.  |

| params setting           | reconstruction    | 
| ---------------------  | ------------------- |
| polyDegree = 0, wendlandRadius = 0.03, Resoluton = 80          |   <img align="center" src="./res/poly0.png" width="300">   |
| polyDegree = 2, wendlandRadius = 0.03, Resoluton = 80 |   <img align="center" src="./res/poly2.png" width="300">   |
| polyDegree = 0, wendlandRadius = 0.03, Resoluton = 50 |   <img align="center" src="./res/resolution50.png" width="300">   |
| polyDegree = 0, wendlandRadius = 0.01, Resoluton = 80 |   <img align="center" src="./res/wendland01.png" width="300">   |

**Please compare the computation time for constraints computation and MLS using brute force and a spatial index. Use hound.off for the comparison and use the same parameters for the timings with and without use of spatial index (please report the parameters that you used).**
| step                    | brute force         |  spatial index      |
| :---------------------: | :-----------------: | :-----------------: |
| constraints             |   1 seconds   |    0.037 seconds             |
| MLS                     |   23 seconds              |    0.26 seconds             |

**Parameters**: hound.off with Resolution = 80, wnedlandRadius = 0.03, polyDegree = 0



### 4 - Theory Question 1

**Prove that the normal of a surface defined by an implicit function is proportional to its gradient.**

Please show your answer in screenshot/photos (or link to a PDF). Make sure your answer is readable. You can add more pages by modifying the table.

[Solution to Theory Questions](./res/Q4.pdf)


### 5 - Theory Question 2

**Compute the closed-form gradient of the MLS approximation.**

Please show your answer in screenshot/photos (or link to a PDF). Make sure your answer is readable. You can add more pages by modifying the table.

[Solution to Theory Questions](./res/Q4.pdf)


### 6 - Normal-based v.s. point-based reconstruction ("hound.off")
| method       | view 01             | view 02            | comments           | 
| :---------:  | ------------------- | ------------------ | ------------------ |
| point-based  |<img align="center" src="./res/Q6_point_v1.png" width="300">| <img align="center"  src="./res/Q6_point_v2.png" width="300"> | xxx |
| normal-based |<img align="center" src="./res/Q6_normal_2.png" width="300">| <img align="center"  src="./res/Q6_normal_1.png" width="300"> | xxx |

**Description of the method**: Instead of blending between constant values, we blend between functions associated with them [Chen Shen et al. 2004](./res/Q6.pdf). In particular, we need to replace the constant value $f_k$ used in the point-based reconstruction with a function in terms of the normal of $p_k$: $S_k(x) = \phi_k + (x - q_k)^T\hat{n}_k$.

**Comparison with point-based reconstruction**: Point-based method can cause undesirable oscillatory behavior as the evaluation point moves away from the surface [Chen Shen et al. 2004](./res/Q6.pdf), while normal-based method does not have this problem. This difference can be seen from the ear of the hound. It can be seen that in the point-based reconstruction, the hound's ear has some sawtooth, but this behaviour does not occur in the normal-based reconstruction.

### 7 - MLS v.s. Screened Poisson Reconstruction v.s. RIMLS

**No implementation required, you can use [Meshlab](https://www.meshlab.net) to perform the comparisons.**

| model names  | MLS          | Possion             | RIMLS               | 
| :---------:  | :---------:  | :-----------------: | :-----------------: |
| cat          |<img align="center" src="./res/Q3_cat.png" width="300">| <img align="center"  src="./res/Q7_Poisson_cat.png" width="300"> |<img align="center"  src="./res/Q7_RIMLS_cat.png" width="300"> |
| luigi        |<img align="center" src="./res/Q3_luigi.png" width="300">| <img align="center"  src="./res/Q7_Poisson_luigi.png" width="300"> |<img align="center"  src="./res/Q7_RIMLS_luigi.png" width="300"> |
| comments        | Both Poisson and RIMLS perform better than MLS with respect to reconstructing the hat of the luigi. | Poisson method produces the most smooth surface reconstruction. The Poisson method considers all the points at once, and is therefore highly resilient to data noise. That's why the result is more smooth. | The RIMLS method is the best at reconstructing details (e.g. eyes and nose of the cat) from relatively sparse data. RIMLS method combines MLS with robust statistics, and therefore is good at preserving small features.  |

### 8 - PCA normals v.s. provided normals (luigi.off)

| model names  | PCA normals             |  Provided normals       | 
| :---------:  | :---------------------: | :---------------------: |
| luigi        | <img align="center"  src="./res/Q5_pca_luigi.png" width="300"> |<img align="center"  src="./res/Q5_provided_luigi.png" width="300">  |
| comments        | The reconstruction using PCA normals failed at the hat and produced some noises around the legs, but succeeded at most of the parts. | The reconstruction using provided normals has less noise around the corners of the mesh. |