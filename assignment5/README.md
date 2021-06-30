# Assignment 5

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations.   

## Required results

## <u>Optimization</u>

##### Gradient Decent method:

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/grad_dec_opt.png?raw=true" />

As we can see, Gradient Decent does find the minimum value of Rosenbrock function (which is 1,1) but it    takes 9750 iterations which is quite a lot.

##### Newton's Method:

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/newton_opt.png?raw=true" />

The Newton's Method does find the minimum also (even more precisely), but as we can see it takes much less iterations (better performance) - as we expected and saw in class.

## <u>Spring system simulation</u>

##### Spring simulation with Gradient Descent method:

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/grad_dec_spring.png?raw=true" />

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/grad_dec_spring_energy.png?raw=true" />

We calculate the opposite of the Force function, at each point of the current Spring.

The calculated value is stored in the gradient vector, at the place of each node.

As we can see, the minimum is granted but, as we expected, it was very slow; 3468 iterations which is quite a lot.

##### Spring simulation with Newton's Method:

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/newton_spring.png?raw=true" />

!<img src="https://github.com/HaifaGraphicsCourses/geometryprocessing2021-sagigvili/blob/master/assignment5/images/newton_spring_energy.png?raw=true" />

Newton's Method gets the same desired result but with much faster converge time - only 3 iterations.