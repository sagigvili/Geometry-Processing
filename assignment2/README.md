# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
**1) Visualization of the constrained points for the 'cat.off' point cloud:**

![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/1.png?raw=true)

**2) Grid with nodes colored according to their implicit function values:**

​		Parameters for both grids: 

​			Resolution - 30

​			polyDegree - 0

​			wendlandRadius - 0.1

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/2.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/3.png?raw=true)

3) Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree:

​		Parameters set #1:

​			Resolution - 20

​			polyDegree - 0

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/4.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/5.png?raw=true)



​		Parameters set #2:

​			Resolution - 30

​			polyDegree - 0

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/6.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/7.png?raw=true)



​		Parameters set #3:

​			Resolution - 40

​			polyDegree - 0

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/8.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/9.png?raw=true)





​		Parameters set #4:

​			Resolution - 20

​			polyDegree - 1

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/10.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/11.png?raw=true)



​		Parameters set #4:

​			Resolution - 20

​			polyDegree - 2

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/12.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/13.png?raw=true)



​		Parameters set #5:

​			Resolution - 20

​			polyDegree - 0

​			wendlandRadius - 0.2	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/14.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/15.png?raw=true)



​		Parameters set #6:

​			Resolution - 20

​			polyDegree - 0

​			wendlandRadius - 0.3	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/16.png?raw=true)

​			

​		Parameters set #7, anisotropic res in the 3 axes :

​			Resolution - 30 x 25 x 20

​			polyDegree - 0

​			wendlandRadius - 0.1	

​		cat.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/17.png?raw=true)



​		luigi.off:

​			![](https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment2/images/18.png?raw=true)

Brief summarize of observations:

As we can see the more we increase resolution, the better reconstructed model we get; when the axes are not homogenous the models getting twitched

Moreover, polyDegree of 1 and 2 increases smoothness of the model and adds a lot of noises around - result of overfitting 

At last, higher wendlandRadius makes the model much more smoother so it most of this parts disappeared

All reconstructed models' off files on results folder		

4) Theory question: Save your notes to assignment2/results and add a link to this page - added as "Theory Question"

### Optional Tasks

1) Added to assignment2/results/Optional Task.docx

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with this method, and try to understand the differences. Report your findings.
