# Assignment 5

Name: 'Di Zhuang'

Legi-Nr: '21958772'

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Tasks

1) **Multiresolution Mesh Editing**: Provide screenshots for 4 different deformed meshes (woody-lo, woody-hi, hand, and cylinder). For each example, provide a rendering of S, B, B' and S'. (questions 1.1 - 1.4 of the assignment sheet)

2) **Real time mesh editing**: Provide animated gifs or short videos for 4 different deformed meshes (bar, bumpy_plane, camel_head, and cactus) showing that your algorithm can perform in real time. (question 1.5 of the assignment sheet)

3) **Deformation transfer**: Discuss and show the differences to the results obtained with the high-frequency detail transfer from part 1.4. on 4 different meshes (bar, bumpy_plane, camel_head, and cactus). (part 2 of the assignment sheet)



## Reports
### 1.1 - 1.4 Multiresolution Mesh Editing
| model name     | S     |  B    |  B'   |  S'   |
| :-----------:  | ----- | ----- | ----- | ----- |
| woody-lo       |<img align="center" src="./res/woody_lo_S.png" width="300">| <img align="center"  src="./res/woody_lo_B.png" width="300"> |<img align="center" src="./res/woody_lo_Bp.png" width="300">| <img align="center"  src="./res/woody_lo_Sp.png" width="300"> |
| woody-hi       |<img align="center" src="./res/woody_hi_S.png" width="300">| <img align="center"  src="./res/woody_hi_B.png" width="300"> |<img align="center" src="./res/woody_hi_Bp.png" width="300">| <img align="center"  src="./res/woody_hi_Sp.png" width="300"> |
| hand           |<img align="center" src="./res/hand_S.png" width="300">| <img align="center"  src="./res/hand_B.png" width="300"> |<img align="center" src="./res/hand_Bp.png" width="300">| <img align="center"  src="./res/hand_Sp.png" width="300"> |
| cylinder       |<img align="center" src="./res/cylinder_S.png" width="300">| <img align="center"  src="./res/cylinder_B.png" width="300"> |<img align="center" src="./res/cylinder_Bp.png" width="300">| <img align="center"  src="./res/cylinder_Sp.png" width="300"> |

### 1.5 Real time mesh editing

Show real time mesh editing using animated gifs or short videos. *Max 15 seconds per gif, better if 5 to 10 seconds*.

| model name     |   S' - real time   |
| :-----------:  | -----  |
| bar            | <img align="center"  src="./res/bar.gif" width="300"> |
| bumpy_plane    | <img align="center"  src="./res/bumpy_plane.gif" width="300"> |
| camel_head     | <img align="center"  src="./res/camel_head.gif" width="300"> |
| cactus         | <img align="center"  src="./res/cactus.gif" width="300"> |


### 2 Deformation transfer
| model name     | High-freq detail transfer             | Deformation transfer                 |
| :-----------:  | ------------------------------------- |------------------------------------- |
| bar            |<img align="center" src="./res/bar_Sp_hf.png" width="300">| <img align="center"  src="./res/bar_Sp_d.png" width="300"> |
| bumpy_plane    |<img align="center" src="./res/bumpy_plane_Sp_hf.png" width="300">| <img align="center"  src="./res/bumpy_plane_Sp_d.png" width="300"> |
| camel_head     |<img align="center" src="./res/camel_head_Sp_hf.png" width="300">| <img align="center"  src="./res/camel_head_Sp_d.png" width="300"> |
| cactus         |<img align="center" src="./res/cactus_Sp_hf.png" width="300">| <img align="center"  src="./res/cactus_Sp_d.png" width="300"> |


#### Observations

|      | High-freq detail transfer             | Deformation transfer                 |
| :-----------:  | ------------------------------------- |------------------------------------- |
| Your Comments  |xxxxxxxxxxxxxxxxxxx                    | It can be seen from the above examples that deformation transfer performs better at preserving the smoothness of the deformed mesh. More importantly, deformation transfer does not lead to self-intersection, while high-frequency detail transfer method may lead to self-intersections when the surface can locally not be properly approximated by a height field.         |

