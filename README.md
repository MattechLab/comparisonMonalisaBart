# Comparison of Monalisa and BART Reconstructions

This repository provides a comparative analysis of Monalisa and BART reconstructions on four different 2D images. The goal is not to establish superiority of any framework but rather to demonstrate that Monalisa achieves results comparable to other well-known open-source frameworks.

We compare both frameworks using L1 and L2 regularization-based reconstructions. To objectively assess the reconstruction quality, we generate synthetic raw MRI measurements, denoted as \( y \), starting from known ground truth images. These synthetic measurements serve as raw data for reconstruction using both frameworks. To ensure fair comparison, we perform a grid search to determine the optimal regularization parameters for each framework. The undersampling strategy follows a 2D radial trajectory with 30 lines of 512 points each, a challenging scenario that highlights the impact of regularization.

## Repository Structure

- `/images/` - Contains the four different test images. The file `radial_2D_CINE.mat` also includes coil sensitivity maps, which remain the same across all images and reconstruction frameworks.
- `generateBARTdatafiles.m` - A script that generates `.cfl` and `.hdr` files in the correct format expected by BART.
- `run_recons_Bart.ipynb` - A jupyter notebook where we run all the bart reconstruction gridsearching for the optimal regularization value.
- `run_recons_Monalisa.m` - A MATLAB script where we run all the Monalisa reconstruction gridsearching for the optimal regularization value.
- `helpers.py` - Bulk of the python code to do the analysis.
- `lineSearchrecon.m` - Bulk of the MATLAB code to do the analysis.
## Reconstruction Methodology

The regularized reconstruction problem is formulated as:

\[
x^* = \arg\min_x \| A x - y \|_2^2 + \lambda R(x)
\]

where the regularization term \( R(x) \) influences the reconstructed image magnitude. Since Structural Similarity Index (SSIM) and L2-distance metrics are sensitive to image scaling, we rescale the reconstructed images post-reconstruction. This rescaling is performed to match the mean intensity of the ground truth image within a manually selected elliptical region of interest (ROI).

### Steps of Comparison

1. **Ground Truth & ROI Selection**: We assume a known coil sensitivity map \( C \) and use four ground truth images: a phantom, a brain slice, and two cardiac slices. An elliptical ROI is manually selected for each image.
2. **Trajectory Generation**: A 2D radial trajectory is generated using:
t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
3. **Data Simulation**: Raw MRI measurements are simulated using:
y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
The corresponding volume elements or density compensation function weights are also computed.
4. **Reconstruction & Evaluation**:
- L1- and L2-regularized iterative reconstructions are performed for both frameworks, using 140 iterations (ensuring convergence).
- The reconstructed images are rescaled to match the mean intensity of the ground truth within the ROI.
- SSIM is computed for various regularization values, and the best-performing parameter (highest SSIM) is selected.

## Results

The best SSIM values for each of the four images are reported for both reconstruction frameworks, providing a quantitative comparison of their performance.
