import scipy.io
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
from scipy.io import loadmat
from scipy.linalg import norm
import json
from datetime import datetime
from bart.python import cfl  # Assuming the provided code is in a module called bart.py
from matplotlib.patches import Ellipse
import ipywidgets as widgets
from IPython.display import display

def readMatToNpArray(matfilepath):
    # Load the .mat file
    mat_data = scipy.io.loadmat(matfilepath)
    # Extract the 'x' variable
    x = mat_data['x']
    # Convert to a NumPy array (if not already)
    x_arr = np.array(x)
    return x_arr

import numpy as np
from numpy.linalg import norm
from skimage.metrics import structural_similarity as ssim

def align_gt_to_recon_magnitude(gt, recon, mask=None):
    """
    Align ground-truth to reconstruction scale 
    using least-squares affine mapping in the *magnitude (real) space*.
    """
    if mask is None:
        mask = np.ones_like(gt, dtype=bool)

    gt_vals = np.abs(gt[mask].ravel())
    recon_vals = np.abs(recon[mask].ravel())

    # Linear regression in real domain
    # recon ≈ a*gt + b
    cov = np.mean((gt_vals - gt_vals.mean()) * (recon_vals - recon_vals.mean()))
    var = np.mean((gt_vals - gt_vals.mean())**2)
    a = cov / var
    b = recon_vals.mean() - a * gt_vals.mean()

    gt_aligned = a * np.abs(gt) + b   # affine transform applied to magnitudes
    return gt_aligned, a, b

def compute_metrics(recon_images, groundtruths, mask=None, percentiles=(0.5, 99.5)):
    """
    Compute SSIM and L2 distance between reconstruction(s) and ground truth(s).
    - Works for a single pair of arrays or a list/iterable of pairs.
    - Aligns ground truth to reconstruction in magnitude space before computing metrics.
    """
    
    # Helper to process one pair
    def _compute_single(recon, gt):
        mag_gt_aligned, _, _ = align_gt_to_recon_magnitude(gt, recon, mask=mask)
        mag_recon = np.abs(recon)

        # Choose ROI for robust dynamic range
        if mask is None:
            roi = mag_gt_aligned
        else:
            roi = mag_gt_aligned[mask]

        lo, hi = np.percentile(roi, percentiles)
        L = hi - lo

        # SSIM
        ssim_val = ssim(mag_recon, mag_gt_aligned, data_range=L)

        # L2 distance
        l2_distance = norm(mag_recon - mag_gt_aligned)

        return ssim_val, l2_distance

    # Case 1: Single image pair
    if isinstance(recon_images, np.ndarray) and isinstance(groundtruths, np.ndarray):
        ssim_val, l2_distance = _compute_single(recon_images, groundtruths)
        return ssim_val, l2_distance

    # Case 2: List or iterable of image pairs
    else:
        ssim_scores, l2_distances = [], []
        for recon, gt in zip(recon_images, groundtruths):
            ssim_val, l2_distance = _compute_single(recon, gt)
            ssim_scores.append(ssim_val)
            l2_distances.append(l2_distance)
        return ssim_scores, l2_distances


# Function to normalize an image based on the ROI
def normalize_within_roi(image, groundtruth, roi_params, debug = False):
    """
    Normalize the image by the mean within the elliptical ROI based on the ground truth.
    
    Parameters:
    - image (ndarray): The reconstructed image.
    - groundtruth (ndarray): The ground truth image.
    - roi_params (dict): Ellipse parameters including center_x, center_y, major_axis, minor_axis, and angle.

    Returns:
    - ndarray: Normalized image.
    """
    image = abs(image)
    groundtruth = abs(groundtruth)
    
    mask = np.zeros_like(groundtruth, dtype=bool)
    y, x = np.ogrid[:groundtruth.shape[0], :groundtruth.shape[1]]
    
    # Ellipse equation: ((x - center_x) * cos(angle) + (y - center_y) * sin(angle))^2 / major_axis^2 
    #                + ((x - center_x) * sin(angle) - (y - center_y) * cos(angle))^2 / minor_axis^2 <= 1
    ellipse = (
        ((x - roi_params["center_x"]) * np.cos(np.deg2rad(roi_params["angle"])) +
            (y - roi_params["center_y"]) * np.sin(np.deg2rad(roi_params["angle"]))) ** 2 / roi_params["major_axis"] ** 2 +
        ((x - roi_params["center_x"]) * np.sin(np.deg2rad(roi_params["angle"])) -
            (y - roi_params["center_y"]) * np.cos(np.deg2rad(roi_params["angle"]))) ** 2 / roi_params["minor_axis"] ** 2 <= 1
    )
    mask[ellipse] = True
    # Debugging: Display the mask (It's correct!)
    if debug:
        plt.figure(figsize=(6, 6))
        plt.imshow(mask, cmap='gray')
        plt.title('ROI Mask')
        plt.axis('off')
        plt.show()
        
    mean_gt = np.mean(groundtruth[mask])
    mean_recon = np.mean(image[mask])
    return image * (mean_gt / mean_recon)



'''
def compute_ssim_magnitude(gt, recon, mask=None, percentiles=(0.5, 99.5)):
    """
    Compute SSIM between reconstruction and GT for MRI data.
    Alignment is performed in magnitude (real space), 
    hence the gt mag is first aligned to the reconstruction mag
    with a linear tranformation, then the SSIM is computed.
    """
    # Align GT magnitude to reconstruction magnitude
    mag_gt_aligned, a, b = align_gt_to_recon_magnitude(gt, recon, mask=mask)

    mag_recon = np.abs(recon)

    if mask is None:
        roi = mag_gt_aligned
    else:
        roi = mag_gt_aligned[mask]

    lo, hi = np.percentile(roi, percentiles)
    L = hi - lo
    
    ssim_val = ssim(mag_recon, mag_gt_aligned, data_range=L) # Data range uses quantiles to be more outlier robust
    return ssim_val
'''
def evaluate_bart_reconstruction_with_roi(regvals, prefix, regtype="l2", ellipse_params=None, iterations = 160):
    """
    Evaluate the performance of BART reconstructions over a range of regularization parameters with ROI-based normalization.
    1. Reconstruction
    2. Normalize the recon to have the same average on the ROI then the ground truth
    3. Compute aligned metrics (we first align the ground truth and then compute SSIM and L2 distance) 

    Parameters:
    - regvals (array-like): Search space for regularization parameters.
    - prefix (str): Prefix for input files and results.
    - regtype (str): Regularization type ("l1" or "l2"). Default is "l2".
    - ellipse_params (dict): Ellipse parameters to define the ROI.

    Returns:
    - dict: Results including metrics and commands, and the selected ROI parameters.
    """
    assert regtype in ["l1", "l2"], "regtype must be either 'l1' or 'l2'."
    
    # Set default ellipse parameters if not provided
    if ellipse_params is None:
        if prefix == 'phantom':
            ellipse_params = {"center_x": 129, "center_y": 127, "major_axis": 88, "minor_axis": 117.5, "angle": 0}
        elif prefix == 'eye':
            ellipse_params = {"center_x": 127, "center_y": 85, "major_axis": 85, "minor_axis": 110, "angle": 0}
        elif prefix == 'cardiac':
            ellipse_params = {"center_x": 143, "center_y": 132.5, "major_axis": 64, "minor_axis": 45, "angle": 68}
        else:   
            raise ValueError(f"Invalid prefix provided: {prefix}")

    # Load the ground truth image: need to adjust the naming convention
    groundtruth = cfl.readcfl(f'./bart_data/image_bart_{prefix}')
    # Initialize arrays to store metrics and commands
    l2_distances = []
    ssim_values = []
    bart_commands = []

    # Loop over the regularization parameters
    for regul_param in regvals:
        if regtype == 'l2':
            bart_command = (
                f"./bart/bart pics -r{regul_param} -l2 -i{iterations} -t ./bart_data/traj_bart_{prefix} "
                f"-p ./bart_data/weights_bart_{prefix}_scaled ./bart_data/kspace_bart_{prefix} ./bart_data/C_bart_{prefix} {regtype}reconweight"
            )
        elif regtype == 'l1':
            bart_command = ( 
                f"./bart/bart pics -R T:3:0:{regul_param} -i{iterations} -m -u {10*regul_param} -C4 -t ./bart_data/traj_bart_{prefix} "
                f"-p ./bart_data/weights_bart_{prefix}_scaled ./bart_data/kspace_bart_{prefix} ./bart_data/C_bart_{prefix} {regtype}reconweight"
            )
#                f"./bart/bart pics -R I:0:{regul_param} -i{iterations} -m -u {10*regul_param} -t ./bart_data/traj_bart_{prefix} "
#                f"-p ./bart_data/weights_bart_{prefix}_scaled ./bart_data/kspace_bart_{prefix} ./bart_data/C_bart_{prefix} {regtype}reconweight"
        print(bart_command)
        subprocess.run(bart_command, shell=True)
        reconstructed = cfl.readcfl(f'{regtype}reconweight')
        reconstructed = normalize_within_roi(reconstructed, groundtruth, ellipse_params) # This is not necessary it's just for display purposes
        plt.imshow(abs(reconstructed), cmap="gray")  # Use grayscale colormap
        plt.colorbar()  # Optional: Show color scale
        plt.axis("off")  # Hide axes
        plt.show()
        ssimval,l2dist = compute_metrics(reconstructed, groundtruth)
        l2_distances.append(l2dist)
        ssim_values.append(ssimval)
        #l2_distances.append(float(np.sqrt(np.sum((abs(groundtruth) - abs(reconstructed)) ** 2))))
        #ssim_values.append(compute_ssim_magnitude(groundtruth,reconstructed)) # This is not simply an SSIM it also contains the linear transformation
        #ssim_values.append(float(ssim(abs(groundtruth), abs(reconstructed), data_range=abs(groundtruth).max() - abs(groundtruth).min())))
        bart_commands.append(bart_command)
    
    bestregval = regvals[np.argmax(ssim_values)]
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"results_{prefix}_{timestamp}.json"
    results = {
        "roi_parameters": ellipse_params,
        "regularization_values": list(regvals),
        "l2_distances": l2_distances,
        "ssim_values": ssim_values,
        "bart_commands": bart_commands,
    }

    #with open(filename, "w") as file:
    #    json.dump(results, file, indent=4)
    #print(f"Results saved to '{filename}'.")
    
    # Plot metrics in separate subplots
    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Plot L2 distance
    axs[0].plot(regvals, l2_distances, '-o', label='L2 Distance', linewidth=1.5, color='blue')
    axs[0].set_xscale('log')
    axs[0].set_ylabel('L2 Distance', fontsize=12)
    axs[0].set_title('L2 Distance vs Regularization Parameter', fontsize=14)
    axs[0].grid(True)
    axs[0].legend(fontsize=10)

    # Plot SSIM
    axs[1].plot(regvals, ssim_values, '-s', label='SSIM', linewidth=1.5, color='red')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Regularization Parameter (-r)', fontsize=12)
    axs[1].set_ylabel('SSIM', fontsize=12)
    axs[1].set_title('SSIM vs Regularization Parameter', fontsize=14)
    axs[1].grid(True)
    axs[1].legend(fontsize=10)

    # Add an overall title
    fig.suptitle(f"Evaluation Metrics Across Regularization Parameters {regtype} {prefix} best regval: {bestregval}", fontsize=16, fontweight='bold')
    
    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

    # Get the corresponding regval
    max_index = np.argmax(ssim_values)
    bestreg = regvals[max_index]
    
    print(f" the best regularization value seems to be: {bestreg}")

    return results



# Global variable to store ellipse parameters
ellipse_parameters = {}

def select_roi(image):
    # Create sliders for the ellipse parameters
    slider_center_x = widgets.IntSlider(value=image.shape[1] // 2, min=0, max=image.shape[1], step=1, description='Center X:')
    slider_center_y = widgets.IntSlider(value=image.shape[0] // 2, min=0, max=image.shape[0], step=1, description='Center Y:')
    slider_major = widgets.IntSlider(value=image.shape[1] // 4, min=1, max=image.shape[1] // 2, step=1, description='Major Axis:')
    slider_minor = widgets.IntSlider(value=image.shape[0] // 4, min=1, max=image.shape[0] // 2, step=1, description='Minor Axis:')
    slider_angle = widgets.IntSlider(value=0, min=0, max=360, step=1, description='Angle (°):')

    # Button to validate
    button_validate = widgets.Button(description="Validate", button_style='success')
    
    # Output area to display the ellipse
    output = widgets.Output()

    # Function to update the plot based on slider values
    def update_plot(cx, cy, major, minor, theta):
        with output:
            output.clear_output(wait=True)
            fig, ax = plt.subplots(figsize=(6, 6))
            
            # Display the image
            ax.imshow(image, cmap='gray')
            
            # Add the ellipse
            ellipse = Ellipse((cx, cy), width=2*major, height=2*minor, angle=theta, edgecolor='red', facecolor='none', lw=2)
            ax.add_patch(ellipse)
            
            ax.set_xlim(0, image.shape[1])
            ax.set_ylim(image.shape[0], 0)
            ax.set_title("Adjust the Ellipse Parameters")
            plt.show()

    # Callback for validate button
    def validate_callback(b):
        # Update the global ellipse parameters dictionary
        ellipse_parameters['center_x'] = slider_center_x.value
        ellipse_parameters['center_y'] = slider_center_y.value
        ellipse_parameters['major_axis'] = slider_major.value
        ellipse_parameters['minor_axis'] = slider_minor.value
        ellipse_parameters['angle'] = slider_angle.value
        
        with output:
            print(f"Ellipse Parameters: Center=({ellipse_parameters['center_x']}, {ellipse_parameters['center_y']}), "
                  f"Major Axis={ellipse_parameters['major_axis']}, Minor Axis={ellipse_parameters['minor_axis']}, "
                  f"Angle={ellipse_parameters['angle']}°")

    # Link sliders to the update function
    widgets.interactive(update_plot,
                        cx=slider_center_x,
                        cy=slider_center_y,
                        major=slider_major,
                        minor=slider_minor,
                        theta=slider_angle)

    # Link button to validation function
    button_validate.on_click(validate_callback)

    # Display the UI
    display(widgets.VBox([slider_center_x, slider_center_y, slider_major, slider_minor, slider_angle, button_validate, output]))

    # Initial plot (ensures the plot is displayed when the function is called)
    update_plot(slider_center_x.value, slider_center_y.value, slider_major.value, slider_minor.value, slider_angle.value)

    # Return the final ellipse parameters once the button is clicked
    return ellipse_parameters