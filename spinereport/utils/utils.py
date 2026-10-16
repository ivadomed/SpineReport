import numpy as np
import platform
from skimage import measure
from scipy.spatial import KDTree
from spinereport.utils.image import Image, zeros_like
from scipy import interpolate, ndimage
from skimage import measure, morphology

def find_symmetry_vector_binary(mask, center, angle_step_deg=1.0, refine_window_deg=2.0, refine_step_deg=0.2):
    """
    Find the symmetry axis of a binary 2D mask constrained to pass through a given center.

    Parameters:
        mask: 2D binary array-like
        center: (row, col) tuple for axis anchor
        angle_step_deg: coarse search step in degrees over [0, 180)
        refine_window_deg: local refinement half-window around coarse optimum
        refine_step_deg: local refinement step in degrees

    Returns:
        angle_rad: axis angle in image coordinates (row, col)
        vector_2d: unit direction vector [drow, dcol]
        min_sum: minimum mismatch sum between pos_coords and mirrored neg_coords
    """
    mask = np.asarray(mask).astype(bool)
    if mask.ndim != 2:
        raise ValueError("mask must be a 2D array")

    coords = np.argwhere(mask)
    if coords.size == 0:
        return 0.0, np.array([1.0, 0.0], dtype=np.float64), 0.0

    h, w = mask.shape
    center = np.asarray(center, dtype=np.float64)
    coords_f = coords.astype(np.float64)
    rel = coords_f - center[None, :]

    def _sum_for_angle(theta):
        drow = np.cos(theta)
        dcol = np.sin(theta)
        normal = np.array([-dcol, drow], dtype=np.float64)
        sym = np.array([drow, dcol], dtype=np.float64)

        signed_dist = rel @ normal
        pos_coords = coords[signed_dist >= 0]
        neg_coords = coords[signed_dist < 0]

        # Project neg_coords in the sym-normal plane
        neg_coords_x = np.round(np.dot(neg_coords, normal)).astype(int)
        neg_coords_y = np.round(np.dot(neg_coords, sym)).astype(int)

        # Project pos_coords in the sym-normal plane
        pos_coords_x = np.round(np.dot(pos_coords, normal)).astype(int)
        pos_coords_y = np.round(np.dot(pos_coords, sym)).astype(int)

        if pos_coords_x.size < 3 or neg_coords_x.size < 3:
            return np.inf

        # Substract min along x
        min_x = np.min(pos_coords_x)
        pos_coords_x -= min_x
        neg_coords_x -= min_x
        
        # Substract min along y
        min_y = np.min(np.concatenate([neg_coords_y, pos_coords_y]))
        pos_coords_y -= min_y
        neg_coords_y -= min_y

        # Inverse neg_coords_x to mirror them
        neg_coords_x = -neg_coords_x

        seg1 = np.zeros((np.max(np.concatenate([pos_coords_x, neg_coords_x]))+1, np.max(np.concatenate([neg_coords_y, pos_coords_y]))+1))
        seg2 = np.zeros((np.max(np.concatenate([pos_coords_x, neg_coords_x]))+1, np.max(np.concatenate([neg_coords_y, pos_coords_y]))+1))
        if pos_coords.size != 0:
            seg1[pos_coords_x, pos_coords_y] = 1
        if neg_coords.size != 0:
            seg2[neg_coords_x, neg_coords_y] = 1
        seg = seg1 - seg2
        diff_sum = np.sum(np.abs(seg))
        return diff_sum

    coarse_angles_deg = np.arange(0.0, 180.0, max(angle_step_deg, 0.1))
    coarse_angles = np.deg2rad(coarse_angles_deg)
    coarse_sums = np.array([_sum_for_angle(theta) for theta in coarse_angles], dtype=np.float64)
    best_idx = int(np.argmin(coarse_sums))
    best_angle = float(coarse_angles[best_idx])
    best_sum = float(coarse_sums[best_idx])

    if refine_step_deg > 0 and refine_window_deg > 0:
        local_deg = np.arange(-refine_window_deg, refine_window_deg + refine_step_deg, refine_step_deg)
        local_angles = best_angle + np.deg2rad(local_deg)
        local_sums = np.array([_sum_for_angle(theta) for theta in local_angles], dtype=np.float64)
        local_best = int(np.argmin(local_sums))
        best_angle = float(local_angles[local_best])
        best_sum = float(local_sums[local_best])

    best_angle = float(np.mod(best_angle, np.pi))
    vector_2d = np.array([np.cos(best_angle), np.sin(best_angle)], dtype=np.float64)
    return best_angle, vector_2d, best_sum

def straighten_coordinates(centerline, spine_centerline, radius, debug_volume=None, debug_tag=''):
    # First pass: per-slice origin, superior, raw anterior (left deferred until anterior is regularized)
    base_list = []
    for i, (origin, normal) in enumerate(zip(centerline['position'].T, centerline['derivative'].T)):
        superior = normal / np.linalg.norm(normal)  # in the S direction

        # Find anterior vector using spine centerline
        anterior_point = intersect_centerline_plane(spine_centerline, origin, superior)
        if anterior_point is None:
            base_list.append(None)
            continue

        anterior = anterior_point - origin
        anterior = anterior / np.linalg.norm(anterior)  # in the A direction

        base_list.append([origin, anterior, superior])

    # # Regularize the anterior vectors across contiguous valid slices to damp noise
    # # from the spine centerline intersection before computing the cross product.
    # anterior_raw = {i: np.array(b[1], copy=True) for i, b in enumerate(base_list) if b is not None}
    # valid_idx = sorted(anterior_raw.keys())
    # if valid_idx:
    #     # Smooth each component with a moving average along the slice axis.
    #     window = 11
    #     kernel = np.ones(window) / window
    #     pad = window // 2
    #     raw = np.array([anterior_raw[i] for i in valid_idx])  # shape (n_valid, 3)
    #     padded = np.pad(raw, ((pad, pad), (0, 0)), mode='edge')
    #     smoothed = np.stack(
    #         [np.convolve(padded[:, k], kernel, mode='valid') for k in range(3)],
    #         axis=1,
    #     )
    #     # Re-project onto the plane perpendicular to the local superior and renormalize,
    #     # so the smoothed anterior stays a valid unit vector of the local base.
    #     for j, i in enumerate(valid_idx):
    #         origin, _, superior = base_list[i]
    #         a = smoothed[j]
    #         a = a - np.dot(a, superior) * superior
    #         n = np.linalg.norm(a)
    #         if n < 1e-8:
    #             continue
    #         base_list[i][1] = a / n

    # Second pass: compute left from the regularized anterior and freeze the tuple
    for i, b in enumerate(base_list):
        if b is None:
            continue
        origin, anterior, superior = b
        left = np.cross(superior, anterior)  # in the L direction
        base_list[i] = (origin, left, anterior, superior)

    # Find longest portion of the centerline with complete base available
    longest_start_idx = None
    longest_length = 0
    current_start_idx = None
    for i, base in enumerate(base_list):
        if base is not None:
            if current_start_idx is None:
                current_start_idx = i
        else:
            if current_start_idx is not None:
                length = i - current_start_idx
                if length > longest_length and current_start_idx <= i-1:
                    longest_length = length
                    longest_start_idx = current_start_idx
                current_start_idx = None
    if current_start_idx is not None:
        length = i - current_start_idx
        if length > longest_length and current_start_idx <= i-1:
            longest_length = length
            longest_start_idx = current_start_idx
    if longest_start_idx is None:
        raise ValueError("Cannot project spine centerline on orthogonal planes of the centerline.")    
    
    coordinates = np.zeros((3, 2*radius+1, 2*radius+1, longest_length), dtype=float)
    centerline_start_idx = centerline['position'].T[longest_start_idx,2].astype(int)
    for i, (origin, left, anterior, superior) in enumerate(base_list[longest_start_idx:longest_start_idx+longest_length]):
        steps = np.linspace(start=-radius, stop=+radius, num=2*radius+1)
        coordinates[:, :, :, i] = origin[:, None, None]
        coordinates[:, :, :, i] += steps[None, :, None] * left[:, None, None]
        coordinates[:, :, :, i] += steps[None, None, :] * anterior[:, None, None]

    # TODO: remove - temporary sagittal view with the orthogonal bases overlaid, so we can
    # eyeball whether (origin, left, anterior, superior) look right along the spine.
    if debug_volume is not None:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        sag_slice = debug_volume[debug_volume.shape[0] // 2, :, :]  # shape (ny, nz)
        used_bases = base_list[longest_start_idx:longest_start_idx + longest_length]
        n_arrows = 60
        stride = max(1, len(used_bases) // n_arrows)
        arrow_len = max(5.0, 0.03 * sag_slice.shape[1])

        fig, ax = plt.subplots(figsize=(12, 7))
        ax.imshow(sag_slice.T, origin='lower', cmap='gray', aspect='auto')
        for base in used_bases[::stride]:
            if base is None:
                continue
            origin, left, anterior, superior = base
            y0, z0 = float(origin[1]), float(origin[2])
            for vec, color in [(left, 'red'), (anterior, 'blue'), (superior, 'green')]:
                dy = vec[1] * arrow_len
                dz = vec[2] * arrow_len
                ax.arrow(
                    y0, z0, dy, dz,
                    color=color, head_width=1.2, head_length=1.8,
                    length_includes_head=True, alpha=0.85, linewidth=0.8,
                )
            ax.plot(y0, z0, '.', color='yellow', markersize=2)

        ax.set_xlabel('A-P (voxels)')
        ax.set_ylabel('I-S (voxels)')
        ax.set_title(f'Sagittal view with straightening bases{" (" + debug_tag + ")" if debug_tag else ""}')
        ax.legend(handles=[
            Line2D([0], [0], color='red', label='left (L)'),
            Line2D([0], [0], color='blue', label='anterior (A)'),
            Line2D([0], [0], color='green', label='superior (S)'),
            Line2D([0], [0], marker='.', color='yellow', linestyle='', label='origin'),
        ], loc='upper right', fontsize=8)
        fig.tight_layout()
        suffix = f'_{debug_tag}' if debug_tag else ''
        fig.savefig(f'straightening_bases_sagittal{suffix}.png', dpi=200)
        plt.close(fig)

    return coordinates, centerline_start_idx

def intersect_centerline_plane(centerline, point, normal):
    for i in range(len(centerline['position'].T)-1):
        coord_i = centerline['position'].T[i]
        coord_i1 = centerline['position'].T[i+1]
        if np.dot(coord_i - point, normal) * np.dot(coord_i1 - point, normal) <= 0:
            break
        elif i == len(centerline['position'].T)-2:
            return None
    while np.dot(coord_i - coord_i1, normal) > 1e-3:
        mid = (coord_i + coord_i1) / 2
        if np.dot(mid - point, normal) * np.dot(coord_i - point, normal) <= 0:
            coord_i1 = mid
        else:
            coord_i = mid
    intersec_point = (coord_i + coord_i1) / 2
    return intersec_point

def _properties2d(mask, dim, radius=2):
    """
    Compute shape properties of a single 2D structure (canal OR spinal cord). Accounts for partial volume information.
    :param mask: 2D input mask in uint8 or float (weighted for partial volume) with a single object.
    :param dim: [px, py]: Physical dimension of the image (in mm). X,Y respectively correspond to AP,RL.
    :param radius: Radius of the cylindrical mask used for diameter computation.
    :return: dict of structure-agnostic metrics, or None if the slice is empty.
    """
    # Check if slice is empty
    if np.sum(mask) == 0:
        print('The slice is empty.')
        return None

    # Patch center (slice center of mass is at the patch center by construction of the straightening)
    x_pos = mask.shape[0]//2
    y_pos = mask.shape[1]//2
    centroid = np.array([x_pos, y_pos])

    # Create vector v (AP direction) and orthogonal w (RL direction)
    v = np.array([0, 1])
    v = v / np.linalg.norm(v)
    w = np.array([1, 0])

    # Compute AP diameter along v
    v_mask = cylindrical_mask(shape=mask.shape, p0=centroid, v=v, radius=radius)
    AP_mask = v_mask*mask
    AP_coords = np.argwhere(AP_mask)
    projections = np.dot(AP_coords, v)
    if projections.any():
        diameter_AP = np.mean([np.max(AP_coords[AP_coords[:,0]==row][:,1]) - np.min(AP_coords[AP_coords[:,0]==row][:,1]) for row in np.unique(AP_coords[:, 0])])*dim[0]
    else:
        diameter_AP = 0

    # Compute RL diameter along w
    w_mask = cylindrical_mask(shape=mask.shape, p0=centroid, v=w, radius=radius)
    RL_mask = w_mask*mask
    RL_coords = np.argwhere(RL_mask)
    projections = np.dot(RL_coords, w)
    if projections.any():
        diameter_RL = np.mean([np.max(RL_coords[RL_coords[:,1]==col][:,0]) - np.min(RL_coords[RL_coords[:,1]==col][:,0]) for col in np.unique(RL_coords[:, 1])])*dim[1]
    else:
        diameter_RL = 0

    # R/L symmetry score
    asymmetry_R_L = 2 * np.sum(mask * np.flip(mask, axis=0)) / (np.sum(mask) + np.sum(np.flip(mask, axis=0)) + 1e-8)

    # Area
    area = np.sum(mask) * dim[0] * dim[1]
    right_area = np.sum(mask[:mask.shape[0]//2+1, :]) * dim[0] * dim[1]
    left_area = np.sum(mask[mask.shape[0]//2:, :]) * dim[0] * dim[1]

    # Eccentricity
    if diameter_AP < diameter_RL:
        eccentricity = np.sqrt(1 - diameter_AP**2/diameter_RL**2) if diameter_RL > 0 else 0
    else:
        eccentricity = -np.sqrt(1 - diameter_RL**2/diameter_AP**2) if diameter_AP > 0 else 0

    # Deal with https://github.com/spinalcordtoolbox/spinalcordtoolbox/issues/2307
    if any(x in platform.platform() for x in ['Darwin-15', 'Darwin-16']):
        solidity = -1
    else:
        solidity = compute_solidity_2d(mask)

    return {
        'area': area,
        'right_area': right_area,
        'left_area': left_area,
        'diameter_AP': diameter_AP,
        'diameter_RL': diameter_RL,
        'centroid': centroid,
        'eccentricity': eccentricity,
        'solidity': solidity,  # convexity measure
        'asymmetry_R_L': asymmetry_R_L,
    }

def cylindrical_mask(shape, p0, v, radius):
    """
    Create a 2D binary mask of a 'cylinder' (thick line) along vector `v` passing through `p0`.
    
    Args:
        shape (tuple): Shape of the 2D image (height, width)
        p0 (np.array): A point [y, x] the vector passes through
        v (np.array): Direction vector [vy, vx]
        radius (float): Cylinder radius (in pixels)

    Returns:
        mask (2D np.array): Binary mask with True inside the cylinder
    """
    h, w = shape
    Y, X = np.meshgrid(np.arange(h), np.arange(w), indexing='ij')
    
    # Shift grid by point
    dy = Y - p0[0]
    dx = X - p0[1]
    
    # Normalize direction vector
    v = v / np.linalg.norm(v)
    
    # Compute perpendicular distance to the line (vector projection method)
    # Distance = ||(point - p0) - ((point - p0) · v) * v||
    dot = dx * v[1] + dy * v[0]
    proj_x = dot * v[1]
    proj_y = dot * v[0]
    perp_x = dx - proj_x
    perp_y = dy - proj_y
    dist = np.sqrt(perp_x**2 + perp_y**2)
    
    # Inside mask if distance < radius
    mask = dist < radius
    return mask

def compute_solidity_2d(mask):
    labeled = measure.label(mask)
    props = measure.regionprops(labeled)[0]
    return props.solidity

def project_point_centerline(centerline, ref_coord):
    dist = np.linalg.norm(centerline - ref_coord, axis=1)
    min_idx = np.argmin(dist)
    return centerline[min_idx]    

def get_centerline(seg, smooth=50):
    '''
    Based on https://github.com/spinalcordtoolbox/spinalcordtoolbox/blob/master/spinalcordtoolbox/centerline/core.py

    Extract centerline from canal segmentation using center of mass and interpolate with bspline
    Expect orientation RPI
    '''
    arr = np.array(np.where(seg.data))
    # Loop across SI axis and average coordinates within duplicate SI values
    sorted_avg = []
    for i_si in np.unique(arr[2]):
        sorted_avg.append(arr[:, arr[2] == i_si].mean(axis=1))
    x_mean, y_mean, z_mean = np.array(sorted_avg).T
    z_ref = np.array(range(z_mean.min().astype(int), z_mean.max().astype(int) + 1))

    # Interpolate centerline
    px, py, pz = seg.dim[4:7]
    x_centerline_fit, x_centerline_deriv = bspline(z_mean, x_mean, z_ref, smooth=smooth, pz=pz)
    y_centerline_fit, y_centerline_deriv = bspline(z_mean, y_mean, z_ref, smooth=smooth, pz=pz)

    # Construct output
    arr_ctl = np.array([x_centerline_fit, y_centerline_fit, z_ref])
    arr_ctl_der = np.array([x_centerline_deriv, y_centerline_deriv, np.ones_like(z_ref)])

    # Create centerline dictionary
    centerline = {"position": arr_ctl, "derivative": arr_ctl_der}
    return centerline

def bspline(x, y, xref, smooth, deg_bspline=3, pz=1):
    """
    Copied from https://github.com/spinalcordtoolbox/spinalcordtoolbox/blob/master/spinalcordtoolbox/centerline/curve_fitting.py
    Bspline interpolation.

    The smoothing factor (s) is calculated based on an empirical formula (made by JCA, based on
    preliminary results) and is a function of pz, density of points and an input smoothing parameter (smooth). The
    formula is adjusted such that the parameter (smooth) produces similar smoothing results than a Hanning window with
    length smooth, as implemented in linear().

    :param x:
    :param y:
    :param xref:
    :param smooth: float: Smoothing factor. 0: no smoothing, 5: moderate smoothing, 50: large smoothing
    :param deg_bspline: int: Degree of spline
    :param pz: float: dimension of pixel along superior-inferior direction (z, assuming RPI orientation)
    :return:
    """
    if len(x) <= deg_bspline:
        deg_bspline -= 2
    density = (float(len(x)) / len(xref)) ** 2
    s = density * smooth * pz / float(3)
    # Then, run bspline interpolation
    tck = interpolate.splrep(x, y, s=s, k=deg_bspline)
    y_fit = interpolate.splev(xref, tck, der=0)
    y_fit_der = interpolate.splev(xref, tck, der=1)
    return y_fit, y_fit_der

def fastest_dilation_edt(mask, radius):
    # Calculate distance from the canal (foreground is 1, background is 0)
    # We want the distance of 0s from the 1s
    distance = ndimage.distance_transform_edt(mask == 0)
    
    # The dilation is where the distance from the original mask is <= radius
    # (Note: we use mask==0, so the distance is 0 inside the canal)
    return distance <= radius