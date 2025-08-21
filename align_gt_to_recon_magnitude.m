function [gt_aligned, a, b] = align_gt_to_recon_magnitude(gt, recon, mask)
    % ALIGN_GT_TO_RECON_MAGNITUDE
    % Align ground truth magnitude to reconstruction magnitude
    % using least-squares affine mapping: recon ≈ a * gt + b.
    %
    % Inputs:
    %   gt    - ground truth complex image
    %   recon - reconstructed complex image
    %   mask  - optional logical mask (same size as gt), default = all true
    %
    % Outputs:
    %   gt_aligned - aligned ground truth (magnitude)
    %   a, b       - affine parameters
    
    if nargin < 3 || isempty(mask)
        mask = true(size(gt));
    end

    gt_vals = abs(gt(mask));
    recon_vals = abs(recon(mask));
    
    % Linear regression: recon ≈ a*gt + b
    covar = mean((gt_vals - mean(gt_vals)) .* (recon_vals - mean(recon_vals)));
    var_gt = mean((gt_vals - mean(gt_vals)).^2);

    a = covar / var_gt;
    b = mean(recon_vals) - a * mean(gt_vals);

    % Apply affine transform to full gt magnitude
    gt_aligned = a * abs(gt) + b;
end