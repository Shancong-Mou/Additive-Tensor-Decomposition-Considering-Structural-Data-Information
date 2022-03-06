# Additive-Tensor-Decomposition-Considering-Structural-Data-Information
These are supplementary materials for the paper 'Additive Tensor Decomposition Considering Structural Data Information'.
It includes three videos for better illustration of the simulation/case study results. The codes for simulations studies can be found in Exp1 and Exp2 folder. You can run main.m files in each folder to reproduice the results.
## 4.	Simulation studies for performance evaluation
### 4.1	Simulation study for Example 1
The data for monitoring the crack growing process include <img src="https://latex.codecogs.com/gif.latex?I_1=30 " /> consecutive measurement images of size 40×40. These images form a tensor M in R^(30×40×40). We simulate M by summing up two tensors <img src="https://latex.codecogs.com/gif.latex?X_1 " /> 
 and <img src="https://latex.codecogs.com/gif.latex?X_2 " />  that represent the background and the crack, respectively. Each <img src="https://latex.codecogs.com/gif.latex?X_1 (i,:,:),i∈[I_1 ], " /> is generated by a 2D smooth Gaussian process representing the background. Most values in the image <img src="https://latex.codecogs.com/gif.latex?X_2 (i,:,:) " />  are zeros, and the non-zero values of <img src="https://latex.codecogs.com/gif.latex?X_2 (i,:,:) " />  gradually grows when the index i increases from 1 to 30, representing the crack growing on the wall. These values are generated from i.i.d. <img src="https://latex.codecogs.com/gif.latex?N(0.1,0.1) " />  random variables to represent the random lighting and shadowing conditions.  The first image in the following video illustrates M, the second images illustrate the corresponding images for actual crack which is a continuous line and the third image illustrate the images for the simulated crack under irregularly illuminated conditions.
We then decompose M into two components <img src="https://latex.codecogs.com/gif.latex?X_1 " />  and <img src="https://latex.codecogs.com/gif.latex?X_2 " />  by solving Problem (2). In the ADMM algorithm, the step size is <img src="https://latex.codecogs.com/gif.latex?\eta=0.01 " /> and the tuning parameters are λ_1,1=λ_1,2=1,λ_2,1=10, and λ_2,2=0.08. The fifth images of Figures 6 (a), (b) illustrate the estimated crack of the 20th and the 30th image using the ATD-based method. It is shown that the ATD-based method captures the growth of the whole crack accurately. 
![Output sample](https://github.com/Sean9511/Additive-Tensor-Decomposition-Considering-Structural-Data-Information/blob/master/Res_sim1.gif)

### 4.2	Simulation study for Example 2
The tensor M in Example 2 is generated to simulate the consecutive measurements taken from a thermal camera in a heated surface monitoring process. It also contains 30 images of size 40×40, and it is generated by summing up three tensors <img src="https://latex.codecogs.com/gif.latex?X_1, X_2, and X_3 " />  of the same size that represents the true background, the static hotspot, and the moving hotspot respectively. Among them, each mode 1 slice of the tensor <img src="https://latex.codecogs.com/gif.latex?X_1 " />   is generated from 
 <img src="https://latex.codecogs.com/gif.latex?X(i,:,:)=U_i T_0+(1-U) T_1 " />
where  <img src="https://latex.codecogs.com/gif.latex?T_0 " /> is a 40×40 matrix representing the heating effect of the heating process, <img src="https://latex.codecogs.com/gif.latex?T_1" /> is a matrix of the same size representing the cooling effect and <img src="https://latex.codecogs.com/gif.latex?U" /> is a <img src="https://latex.codecogs.com/gif.latex?U[0,1]" /> random variable representing a random combination of the two effects. 
The images representing the matrices <img src="https://latex.codecogs.com/gif.latex?T_0" /> and <img src="https://latex.codecogs.com/gif.latex?T_1" /> are shown in the following figure. 
![alt text](https://github.com/Sean9511/Additive-Tensor-Decomposition-Considering-Structural-Data-Information/blob/master/Heat%26Cooling.png?raw=true)
To simulate the heating effect of a single point heating source at the center of this image, we generate <img src="https://latex.codecogs.com/gif.latex?T_0 (i,j)" /> using the value of <img src="https://latex.codecogs.com/gif.latex?f_0 (i,j)" />, where <img src="https://latex.codecogs.com/gif.latex?f_0" /> is the probability density function of <img src="https://latex.codecogs.com/gif.latex?N((20,20)^T,10I)" />, where I is a 2×2 identity matrix. Then, we transform all values of <img src="https://latex.codecogs.com/gif.latex?T_0" /> linearly such that the maximum and minimum value of <img src="https://latex.codecogs.com/gif.latex?T_0" /> are 1 and 0, respectively. With this setup, the maximum value within <img src="https://latex.codecogs.com/gif.latex?T_0" /> is 1, located at the center of the image; when the pixel moves farther way from the center, the value of <img src="https://latex.codecogs.com/gif.latex?T_0 (i,j)" /> gradually drops to 0.  To simulate the cooling effect, we generate <img src="https://latex.codecogs.com/gif.latex?T_1 (i,j))" /> using a linear function <img src="https://latex.codecogs.com/gif.latex?f_1 (i,j)=c_1 (i+j)+c_2" />, where <img src="https://latex.codecogs.com/gif.latex?c_1" /> and <img src="https://latex.codecogs.com/gif.latex?c_2" /> are adjusted so that the maximum and minimum value within the matrix <img src="https://latex.codecogs.com/gif.latex?T_1" /> are 0 and 1 respectively. It represents that the coolant for the surface flows from the upper-left corner to the bottom-right corner of the image. 

Each image within the tensor <img src="https://latex.codecogs.com/gif.latex?X_2" /> are the same, and the non-zero values in these images are located in a fixed 2×2 rectangle with intensity value 1 in their lower-left corners. The non-zero values in each image of the tensor <img src="https://latex.codecogs.com/gif.latex?X_3" /> are also located in a 2×2 rectangle with intensity value 1. However, this rectangle locates on the upper-left part of the image.  When the image index i increases, the rectangle in <img src="https://latex.codecogs.com/gif.latex?X_3 (i,:,:)" /> moves from the left side to the right side across the images.

We decompose the tensor M into components <img src="https://latex.codecogs.com/gif.latex?X_1" />, <img src="https://latex.codecogs.com/gif.latex?X_2" /> and <img src="https://latex.codecogs.com/gif.latex?X_3" />  by solving Problem (3). The step size <img src="https://latex.codecogs.com/gif.latex?\eta=0.01" />. Tuning parameters are λ_1,1=λ_1,2=30, λ_1,3= λ_2,1=1, λ_2,2=1.9 and λ_3,1=2.

In the following video, we illustrate the simulated images, with the decomposed background, the static hotspot, and the moving hotspot using the ATD method. The decomposed background and the moving object from the RPCA method are illustrated in the third column. The video illustrating the result of the decomposition is provided in the supplementary material. We can see that the proposed ATD method successfully separates the static background, static hotspot and the moving hotspot. Comparing the result of the ATD, the RPCA can only identify the moving object but fails to separate the true background and the static hotspot. This is because the RPCA does not consider the spatial smoothness of the background. 
![Output sample](https://github.com/Sean9511/Additive-Tensor-Decomposition-Considering-Structural-Data-Information/blob/master/Res_sim2.gif)

## 5.	A case study with the 3D medical image analysis
With tuning parameters λ_1,1=30, λ_2,1=0.7, λ_2,2=0.16 and λ_3,1=1, we solved this problem using the ADMM algorithm. The images in first row of the following video illustrates the non-zero values from the 3rdth-7th images in the solution of <img src="https://latex.codecogs.com/gif.latex?X_2" />, reflecting the AVC regions, that determines the boundary of the 3D-printed prototype of patients’ aortic root anatomies. It shows that the smooth change of the AVC regions across images is captured and it reflects the anatomic reality. The following video illustrates the result of the decomposition.
For comparison, we also used SSD to extract the AVC regions from individual images (the second row of the video). It achieves similar performance in the 3rd-6th images. In the 7th image, however, the SSD method failed to fully extract the anomaly region. The reason is that the intensity of the AVC regions is not so high. Unlike the ATD formulation, SSD does not consider the similarity of the AVC regions across images so that  the pixels with low intensity on the AVC regions cannot be preserved. This case study further illustrates the versatility of the ATD framework and demonstrates its ability to solving real-world problems. 
![Output sample](https://github.com/Sean9511/Additive-Tensor-Decomposition-Considering-Structural-Data-Information/blob/master/Res_case.gif)

