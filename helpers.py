import numpy as np
import subprocess
import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
import json
from datetime import datetime
from bart.python import cfl  # Assuming the provided code is in a module called bart.py
from matplotlib.patches import Ellipse
import ipywidgets as widgets
from IPython.display import display

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

def evaluate_bart_reconstruction_with_roi(regvals, prefix, regtype="l2", ellipse_params=None, iterations = 30):
    """
    Evaluate the performance of BART reconstructions over a range of regularization parameters with ROI-based normalization.

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
            ellipse_params = {"center_x": 258, "center_y": 254, "major_axis": 176, "minor_axis": 235, "angle": 0}
        elif prefix == 'brain':
            ellipse_params = {"center_x": 271, "center_y": 250, "major_axis": 202, "minor_axis": 220, "angle": 0}
        elif prefix == 'cardiac':
            ellipse_params = {"center_x": 250, "center_y": 234, "major_axis": 107, "minor_axis": 95, "angle": 75}
        elif prefix == 'cardiac2':
            ellipse_params = {"center_x": 262, "center_y": 217, "major_axis": 100, "minor_axis": 132, "angle": 145}
        elif prefix == 'cardiac3':
            ellipse_params = {"center_x": 273, "center_y": 268, "major_axis": 128, "minor_axis": 128, "angle": 0}
        else:   
            raise ValueError(f"Invalid prefix provided: {prefix}")

    # Load the ground truth image: need to adjust the naming convention
    groundtruth = cfl.readcfl(f'./bart_data/final{prefix}image')

    # Initialize arrays to store metrics and commands
    l2_distances = []
    ssim_values = []
    bart_commands = []

    # Loop over the regularization parameters
    for regul_param in regvals:
        if regtype == 'l2':
            bart_command = (
                f"./bart/bart pics -r{regul_param} -l2 -i{iterations} -t ./bart_data/traj_bart "
                f"-p ./bart_data/weights_bart_scaled ./bart_data/kspace_bart{prefix} ./bart_data/C_bart {regtype}reconweight"
            )
        elif regtype == 'l1':
            bart_command = (
                f"./bart/bart pics -R I:1:{regul_param} -i{iterations} -m -u {10*regul_param} -t ./bart_data/traj_bart "
                f"-p ./bart_data/weights_bart_scaled ./bart_data/kspace_bart{prefix} ./bart_data/C_bart {regtype}reconweight"
            )
            
        subprocess.run(bart_command, shell=True)
        reconstructed = cfl.readcfl(f'{regtype}reconweight')
        plt.imshow(abs(reconstructed), cmap="gray")  # Use grayscale colormap
        plt.colorbar()  # Optional: Show color scale
        plt.axis("off")  # Hide axes
        plt.show()
        reconstructed = normalize_within_roi(reconstructed, groundtruth, ellipse_params)
        l2_distances.append(float(np.sqrt(np.sum((abs(groundtruth) - abs(reconstructed)) ** 2))))
        ssim_values.append(float(ssim(abs(groundtruth), abs(reconstructed), data_range=abs(groundtruth).max() - abs(groundtruth).min())))
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





