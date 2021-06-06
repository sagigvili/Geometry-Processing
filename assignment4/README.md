# Assignment 4

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment4/results' and refer to or directly show them in this page.

## Required results for assignment 4

### Mandatory Tasks

Provide screenshots for 4 different deformed meshes. For each example, provide a rendering of S, B, B' and S':

Last column in each line represents the efficiency of the Cholesky factorization, that is average time to calculate each frame in the given mesh.

|    Mesh    |                              S                               |                              B                               |                              B'                              |                              S'                              | Average execution time (seconds) |
| :--------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :------------------------------: |
|  woody-lo  | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/woody-lo-S.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/woody-lo-B.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/woody-lo-B-tag.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/woody-lo-S-tag.PNG?raw=true" style="zoom:50%;" /> |              0.023               |
|    hand    | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/hand-S.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/hand-B.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/hand-B-tag.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/hand-S-tag.PNG?raw=true" style="zoom:50%;" /> |                2                 |
|    bar     | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/bar-S.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/bar-B.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/bar-B-tag.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/bar-S-tag.PNG?raw=true" style="zoom:50%;" /> |               0.9                |
| camel_head | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/camel-S.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/camel-B.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/camel-B-tag.PNG?raw=true" style="zoom:50%;" /> | <img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment4/results/camel-S-tag.PNG?raw=true" style="zoom:50%;" /> |               1.6                |

#### Notes

Dragging the handles (like woody-lo hands) far away from its original position, local self intersection on meshes can be observed.

Similar phenomenon happens when dragging with low amount of handles.

### Optional Task

Discuss and show the differences to the results from the previous task. 
