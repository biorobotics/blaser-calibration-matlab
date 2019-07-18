# blaser_calibration_matlab
## Capturing images and data for new calibration process
- Requires `blaser_ros` and `ur5e_driver`

First, start the ur5e driver and blaser processes.
### Bring up UR5 arm
[ur5e bringup info](https://github.com/biorobotics/ur5e_driver/blob/blaser/README.md). You only need to run the `roslaunnch ur5e_bringup.launch`

[blaser bringup info](https://github.com/biorobotics/blaser-ros/blob/blaser_rewrite/README.md).

Under `blaser-ros/blaser_rewrite` use the `add_im.py` script to capture calibration images. Usage:

`python add_im.py /path/to/calibration/directory`

Only 20-30 images are needed for calibration, but more are welcome. (It may slow the calibration significantly though)

## Running calibration process
In the MATLAB script `unified_proc.m` change `n_val` to the number of calibration images captured. Also add your calibration directory to the MATLAB path.

The script should just run from then on. It's highly recommended to run it section by section (AHEM SARA)

The sections are as follows:
 1. Load calibration data files
 2. Initialize camera intrinsics
 3. Optimize camera intrinsics by minimizing reprojection error
 4. Show reprojection errors of checkerboard
 5. Initialize hand-eye
 6. Optimize handeye simultaneously with intrinsics
 7. Show hand-eye fit
 8. Initialize laser plane params
 9. Optimize laser plane and handeye and intrinsics all simultaneously
 10. Show laser error
 11. Display results in a nice way (ready for copy/pasting into the calibration yaml)
