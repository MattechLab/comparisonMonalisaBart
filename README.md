# Comparison of Monalisa and BART Reconstructions

This repository provides a comparative analysis of Monalisa and BART reconstructions on four different 2D images. The goal is not to establish superiority of any framework but rather to demonstrate that Monalisa achieves results comparable to other well-known open-source frameworks.

We compare both frameworks using L1 and L2 regularization-based reconstructions. To objectively assess the reconstruction quality, we generate synthetic raw MRI measurements, denoted as \( y \), starting from known ground truth images. These synthetic measurements serve as raw data for reconstruction using both frameworks. To ensure fair comparison, we perform a grid search to determine the optimal regularization parameters for each framework. The undersampling strategy follows a 2D radial trajectory with 30 lines of 512 points each, a challenging scenario that highlights the positive impact of regularization.

## Repository Structure

- `/images/` - Contains the four different test images. The file `radial_2D_CINE.mat` also includes coil sensitivity maps, which remain the same across all images and reconstruction frameworks.
- `generateBARTdatafiles.m` - A script that generates `.cfl` and `.hdr` files in the correct format expected by BART.
- `run_recons_Bart.ipynb` - A jupyter notebook where we run all the bart reconstruction gridsearching for the optimal regularization value.
- `run_recons_Monalisa.m` - A MATLAB script where we run all the Monalisa reconstruction gridsearching for the optimal regularization value.
- `helpers.py` - Bulk of the python code to do the analysis.
- `lineSearchrecon.m` - Bulk of the MATLAB code to do the analysis.
## Reconstruction Methodology

The regularized reconstruction problem is formulated as:

$$
x^* = \arg\min_x \| A x - y \|_2^2 + \lambda R(x)
$$

where the regularization term \( R(x) \) influences the reconstructed image magnitude. Since Structural Similarity Index (SSIM) and L2-distance metrics are sensitive to image scaling, we rescale the reconstructed images post-reconstruction. This rescaling is performed to match the mean intensity of the ground truth image within a manually selected elliptical region of interest (ROI).

### Steps of Comparison

1. **Ground Truth & ROI Selection**: We assume a known coil sensitivity map \( C \) and use four ground truth images: a phantom, a brain slice, and two cardiac slices. An elliptical ROI is manually selected for each image.
2. **Trajectory Generation**: A 2D radial trajectory is generated using:
t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
3. **Data Simulation**: Raw MRI measurements are simulated using:
y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
The corresponding volume elements or density compensation function weights are also computed.
4. **Reconstruction & Evaluation**:
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

## Results

The best SSIM values for each of the four images are reported for both reconstruction frameworks, providing a quantitative comparison of their performance. 

Images reconstructed with l1 and l2 regularization are presented below respectively. The quantitative metrics are reported accordingly. l2 distance computes the sum of squared differences between corresponding pixel intensities in two images. It is a straightforward, pixel-wise error measure of the reconstruction that treats every pixel equally. SSIM, on the other hand, is designed to mimic the human visual system by considering perceptual phenomena. It evaluates image quality based on three key components: luminance, contrast, and structure. Visually, for cardiac images both frameworks yield similar image quality, while for the brain and the phantom image, BART's reconstructions seems more affected by artifacts. Moreover, the reported metrics reveal a trade-off: when one framework achieves a higher \gls{ssim} score, the other tends to have a lower l2 distance, indicating that neither framework consistently produces superior image reconstructions. 

![l1images](https://github.com/user-attachments/assets/359fd59a-159d-405c-997a-b0c417d07e53)
![l2images](https://github.com/user-attachments/assets/47c91043-3c43-4c18-893c-7a0239ab9bfd)

## Comparison of l1-regularized reconstructions

| Metric          | Gridded Recon | Monalisa l1-Reg | BART l1-Reg |
|----------------|--------------|----------------|-------------|
| **Phantom**    |              |                |             |
| **SSIM**       | 0.3789       | 0.5867         | **0.6883**  |
| **l2-Distance**| 93.1373      | **29.3729**    | 52.0827     |
| **Brain**      |              |                |             |
| **SSIM**       | 0.2391       | 0.5332         | **0.5645**  |
| **l2-Distance**| 28553.3950   | **13523.8670** | 13531.1620  |
| **Cardiac (Close)** |        |                |             |
| **SSIM**       | 0.7227       | 0.8342         | **0.8448**  |
| **l2-Distance**| 28.4789      | **17.9613**    | 17.9725     |
| **Cardiac (Far)** |          |                |             |
| **SSIM**       | 0.4818       | **0.7930**     | 0.7503      |
| **l2-Distance**| 103.8571     | **43.4874**    | 59.4864     |

*Comparison of l1-regularized reconstructions between BART and Monalisa, including Gridded Reconstructions as a baseline. The best values are in **bold**.*

---

## Comparison of l2-regularized reconstructions

| Metric          | Gridded Recon | Monalisa l2-Reg | BART l2-Reg |
|----------------|--------------|----------------|-------------|
| **Phantom**    |              |                |             |
| **SSIM**       | 0.3789       | **0.5049**     | 0.4670      |
| **l2-Distance**| 93.1373      | 47.4290        | **47.0134** |
| **Brain**      |              |                |             |
| **SSIM**       | 0.2391       | 0.5370         | **0.5567**  |
| **l2-Distance**| 28553.3950   | **13309.8610** | 15796.3090  |
| **Cardiac (Close)** |        |                |             |
| **SSIM**       | 0.7227       | 0.8388         | **0.8501**  |
| **l2-Distance**| 28.4789      | 18.0503        | **17.3688** |
| **Cardiac (Far)** |          |                |             |
| **SSIM**       | 0.4818       | **0.8008**     | 0.7438      |
| **l2-Distance**| 103.8571     | **49.7409**    | 60.3805     |

*Comparison of L2-regularized reconstructions between BART and Monalisa, including Gridded Reconstructions as a baseline. The best values are in **bold**.*
