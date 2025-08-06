# Comparison of Monalisa and BART Reconstructions

This repository provides a comparative analysis of Monalisa and BART reconstructions on four different 2D images. The goal is not to establish superiority of any framework but rather to confirm that Monalisa achieves results comparable to BART.

We compare both frameworks using L1 and L2 regularization-based reconstructions. To assess the reconstruction quality, we generate synthetic raw MRI measurements, denoted as \( y \), starting from known ground truth 2D images and coil sensitivities. These synthetic measurements serve as raw data for reconstruction using both frameworks. With the aim of fair comparison, we perform a grid search to determine the optimal regularization parameters for each framework separately. The undersampling strategy follows a 2D radial trajectory with 30 lines of 512 points each, a challenging scenario that highlights the positive impact of regularization.

## Reconstruction Methodology

The regularized reconstruction problem is formulated as:

$$
x^* = \arg\min_x \| A x - y \|_2^2 + \lambda R(x)
$$

where the regularization term \( R(x) \) influences both the obtained image \( x^* \) and the reconstructed image magnitude. Since Structural Similarity Index (SSIM) and L2-distance metrics are sensitive to image scaling, we rescale the reconstructed images post-reconstruction. This rescaling is performed to match the mean intensity of the ground truth image within a manually selected elliptical region of interest (ROI).

### Steps of Comparison

1. **Ground Truth & ROI Selection**: We assume a known coil sensitivity map \( C \) and use three ground truth images: a phantom, an eye image, and a cardiac slice. An elliptical ROI is manually selected for each image.
2. **Trajectory Generation**: A 2D radial trajectory, constructed by rotating each line sequentially by an angle of \(\frac{2\pi}{30}\) radians, is generated using:
t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
The corresponding volume elements or density compensation function weights are also computed using:
ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');

3. **Data Simulation**: Raw MRI measurements are simulated using:
y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);

4. **Data conversion**: The data is inherently in following the convention of Monalisa. To run BART reconstruction we need to convert the data into the BART convention (eg: change volume element definitions and generate correctly formatted  .clf and .hdr files). This is done in the generateBARTfiles.m script.

5. **Reconstruction & Evaluation**:
- L1-regularized
  
$$
x^* = \arg\min_x \| A x - y \|_2^2 + \lambda \| x \|_1
$$
  
  and L2-regularized
  
$$
x^* = \arg\min_x \| A x - y \|_2^2 + \lambda \| x \|_2
$$

  iterative reconstructions are performed for both frameworks, using 140 iterations (ensuring convergence). \lambda is the regularization parameter
- The reconstructed images are rescaled to match the mean intensity of the ground truth within the ROI.
- SSIM is computed for various regularization values, and the best-performing parameter (highest SSIM) is selected and reported.

## Repository Structure

- `/images/` - Contains the three different test images: A cardiac image, a brain image with the FoV centered on the eye, and the well known Shepp-Logan phanom. The files are .mat files containing also the corresponding coil sensitivity maps used for the multicoil simulation. C_simu.mat is a simulated coil sensitivity map, that is smooth, and is used for the phantom data.
- `generateBARTdatafiles.m` - A script that generates `.cfl` and `.hdr` files in the correct format expected by BART.
- `run_recons_Bart.ipynb` - A jupyter notebook where we run all the bart reconstruction gridsearching for the optimal regularization value.
- 
- `run_recons_Monalisa.m` - A MATLAB script where we run all the Monalisa reconstruction gridsearching for the optimal regularization value. Additionally we run the final recontruction with the optimal regularization value and save the results.
- `helpers.py` - Bulk of the python code to do the analysis.
- `lineSearchrecon.m` - Bulk of the MATLAB code to do the analysis.
- `/reconstructions/...` - All the resulting reconsturctions files
- `/results/...` - Nice plots of the resulting reconstructions
- `FinalReconsEvaluation.ipynb` - Runs the final reconstruction for BART and then generates nice plots of the reconstructions comparing both frameworks.

## Results

The best SSIM values for each of the four images are reported for both reconstruction frameworks, providing a quantitative comparison of their performance. 

Images reconstructed with l1 and l2 regularization are presented below respectively. The quantitative metrics are reported accordingly. l2 distance computes the sum of squared differences between corresponding pixel intensities in two images. It is a straightforward, pixel-wise error measure of the reconstruction that treats every pixel equally. SSIM, on the other hand, is designed to mimic the human visual system by considering perceptual phenomena. It evaluates image quality based on three key components: luminance, contrast, and structure. Visually, for eye images both frameworks yield similar image quality, while for the cardiac and especially the phantom l1 regularized image, BART's reconstructions seems more affected by artifacts. Moreover, in the l2 regularized case the reported metrics reveal a trade-off: when one framework achieves a higher \gls{ssim} score, the other tends to have a lower l2 distance, indicating that neither framework consistently produces superior image reconstructions. 



## Comparison of l1-regularized reconstructions

| Metric          | Gridded Recon | Monalisa l1-Reg | BART l1-Reg |
|----------------|--------------|----------------|-------------|
| **Phantom**    |              |                |             |
| **SSIM**       | 0.2492       | **0.7618**     | 0.5152      |
| **l2-Distance**| 31.6767      | **14.1150**    | 35.5927     |
| **Brain**      |              |                |             |
| **SSIM**       | 0.3240       | **0.7261**     | 0.7172      |
| **l2-Distance**| 94.8667      | **36.5786**    | 37.8975     |
| **Cardiac (Close)** |        |                |             |
| **SSIM**       | 0.4882       | **0.8201**     | 0.7458      |
| **l2-Distance**| 31.3764      | **14.7013**    | 20.5866     |

*Comparison of l1-regularized reconstructions between BART and Monalisa, including Gridded Reconstructions as a baseline. The best values are in **bold**.*

---

## Comparison of l2-regularized reconstructions

| Metric          | Gridded Recon | Monalisa l2-Reg | BART l2-Reg |
|----------------|--------------|----------------|-------------|
| **Phantom**    |              |                |             |
| **SSIM**       | 0.2492       | **0.4453**     | 0.3750      |
| **l2-Distance**| 31.6767      | 26.6510        | **26.3846** |
| **Brain**      |              |                |             |
| **SSIM**       | 0.3240       | **0.7239**     | 0.6720      |
| **l2-Distance**| 94.8667      | 42.0751        | **39.6113** |
| **Cardiac (Close)** |        |                |             |
| **SSIM**       | 0.4882       | **0.8158**     | 0.7275      |
| **l2-Distance**| 31.3764      | **15.6514**    | 20.3854     |

*Comparison of L2-regularized reconstructions between BART and Monalisa, including Gridded Reconstructions as a baseline. The best values are in **bold**.*
