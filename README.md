# SpineReport
[![preprint](https://img.shields.io/badge/preprint-arXiv:2606.10021-orange)](https://arxiv.org/abs/2606.10021)
[![Docs](https://img.shields.io/badge/docs-GitHub_Pages-blue)](https://ivadomed.github.io/SpineReport/)

Automated extraction of spinal morphometrics and generation of structured radiological reports from MRI data.

| Report description |
| :---: |
| <img src="spinereport/resources/imgs/report_first_page.png" width="900"> | 

| Canal report |
| :---: |
| <img src="spinereport/resources/imgs/report_canal.png" width="900"> | 

| Discs report |
| :---: |
| <img src="spinereport/resources/imgs/report_discs.png" width="900"> |

| Foramens report |
| :---: |
| <img src="spinereport/resources/imgs/report_foramens.png" width="900"> |

| Vertebrae report |
| :---: |
| <img src="spinereport/resources/imgs/report_vertebrae.png" width="900"> |

| Cerebro Spinal Fluid (CSF) report |
| :---: |
| <img src="spinereport/resources/imgs/report_csf.png" width="900"> |

In these report examples:
- The background violin plot in gray (or lineplot for the canal) corresponds to the **control group** (multiple subjects)
- The red line in each graph and the pictures correspond to a specific subject in the **test group**. Reports are only generated for the test group.

You can find an example report [here](spinereport/resources/imgs/example_report.pdf)

## How to install ?

1. Open a `bash` terminal in the directory where you want to work.

2. Create and activate a virtual environment using python >=3.10 (highly recommended):
   - venv
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
   - conda env
   ```
   conda create -n myenv python=3.10
   conda activate myenv
   ```

3. Install TotalSpineSeg:
   see information [here](https://github.com/neuropoly/totalspineseg#installation)

4. Install SpineReport:
   ```bash
   git clone git@github.com:ivadomed/SpineReport.git
   cd SpineReport
   pip install .
   ```

## How to generate the reports ?

1. Regroup all niftii files into a folder
> They need to follow the BIDS naming convention (i.e. `sub-<participant_label>_<contrast>.nii.gz`)
```
INPUT_FOLDER/
    subject1_T1w.nii.gz
    subject2_T2w.nii.gz
    subject3_T2w.nii.gz
    subject7_T2w.nii.gz
    ...
``` 

2. Run TotalSpineSeg to get generate the segmentations in the 1mm isotropic space
> Perform this for both your control group and test group
```
totalspineseg INPUT_FOLDER TOTALSPINESEG_FOLDER --iso
```

3. Run SpineReport to generate the reports
```
spinereport -t TEST_TOTALSPINESEG_FOLDER -c CONTROL_TOTALSPINESEG_FOLDER -o reports
```

The **test subjects** corresponds to the subjects for which a report will be generated. The **control subjects** corresponds to the violin plots shown in gray in the background (see report examples at the top of the README.). If you want to generate reports for all your subjects/scans, you can specify the same path for both **TEST_TOTALSPINESEG_FOLDER** and **CONTROL_TOTALSPINESEG_FOLDER**.

## Using your own segmentations

SpineReport can also work with segmentations you produced yourself or with other automatic tools (spinal cord, canal, vertebrae, discs), and fall back on TotalSpineSeg for anything you don't provide. Any per-structure folder passed through the flags below overrides its counterpart in `--test-tss` / `--control-tss`; you can mix and match freely. Be careful, the segmentations must be in the same 1mm isotropic space as the input images in TOTALSPINESEG_FOLDER/input. To ensure that all segmentations are in the correct space, you should run TotalSpineSeg with the `--iso` flag on your raw images, then compute your segmentations onto the `TOTALSPINESEG_FOLDER/input` folder and save them in a separate segmentation folder.

You can also skip TotalSpineSeg entirely by providing every field yourself (`--*-images-dir`, `--*-labels-dir`, `--*-sc-seg-dir`, `--*-canal-seg-dir`, `--*-vertebrae-seg-dir`, `--*-discs-seg-dir`), as long as all images and segmentations are resampled to the same 1mm isotropic space.

### Folder layout

All segmentation folders must be **flat** (no per-subject subdirectories) and every file must share the same BIDS-style basename as the raw image (e.g. `<canal-seg-folder>/sub-001_T2w.nii.gz`).

### Segmentation expectations

- **Spinal cord (`--*-sc-seg-dir`)** and **canal (`--*-canal-seg-dir`)** must be **binary** (any non-zero voxel is treated as foreground). The CSF region is derived as *canal minus SC*.
- **Vertebrae (`--*-vertebrae-seg-dir`)** and **discs (`--*-discs-seg-dir`)** must be multi-label (one integer per anatomical level). Each folder must also contain a `map.json` that names every label in your segmentation using the [tss_map.json](https://github.com/neuropoly/totalspineseg/blob/main/totalspineseg/resources/labels_maps/tss_map.json) keys:
```json
// VERTEBRAE/map.json — keys are tss_map.json names, values are YOUR integer labels
{
  "C1": 1, "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6, "C7": 7,
  "T1": 8, "T2": 9, "T12": 19,
  "L1": 20, "L5": 24,
  "sacrum": 25
}
```
```json
// DISCS/map.json
{
  "C2-C3": 1, "C3-C4": 2, "C6-C7": 5,
  "T1-T2": 7, "T12-L1": 18,
  "L1-L2": 19, "L5-S": 23
}
```
SpineReport uses these `map.json` files to remap your segmentation values to the canonical TotalSpineSeg label scheme before extracting metrics.

- **Landmark labels (`--*-labels-dir`)** are single-voxel labels at the posterior tip of each disc (same idea as TotalSpineSeg's `step1_levels/`). Their integer values must follow the [levels_maps.json](https://github.com/neuropoly/totalspineseg/blob/main/totalspineseg/resources/labels_maps/levels_maps.json) convention (`C1=1`, `C1-C2=2`, ..., `L5-S=25`) see [convention](https://spinalcordtoolbox.com/stable/user_section/tutorials/vertebral-labeling/labeling-conventions.html).

### Example

Use custom spinal cord segmentations from [SCT](https://spinalcordtoolbox.com/):

1. First run TotalSpineSeg to obtain the rest of the segmentations using the flag `--iso` (canal, vertebrae, discs, landmarks):

```bash
totalspineseg raw out-iso --iso
```

2. Then run `sct_deepseg` to get spinal cord segmentations for all your images:
> Because all images must be resampled to 1mm isotropic before using SpineReport, we can use the `out-iso/input` folder as input for `sct_deepseg`

```bash
mkdir -p sc-seg
for file in out-iso/input/*; do 
    base=$(basename "$file"); 
    sct_deepseg spinalcord -i "$file" -o sc-seg/${base/_0000/};
done
```

3. Finally run SpineReport by specifying your new spinal cord segmentations

```bash
spinereport -t out-iso -c out-iso --test-sc-seg-dir sc-seg --control-sc-seg-dir sc-seg -o reports
```

## How to generate group analysis ?

It is also possible to only generate group analysis (by sex and age) with this repository

1. Regroup all niftii files into a folder
> They need to follow the BIDS naming convention (i.e. `sub-<participant_label>_<contrast>.nii.gz`)
```
INPUT_FOLDER/
    subject1_T1w.nii.gz
    subject2_T2w.nii.gz
    subject3_T2w.nii.gz
    subject7_T2w.nii.gz
    ...
```

2. Construct a `tsv` file with demographics like this DEMOGRAPHICS.tsv
> The order of the columns is not important

| participant_id | sex | age |
| :---: | :---: | :---: |
| sub-001 | M | 39 |
| sub-004 | F | 25 |
| sub-006 | F | 34 |

3. Run this command
```
spinereport_plot_by_group -i TOTALSPINESEG_FOLDER -d DEMOGRAPHICS.tsv -o group_analysis
```

### Examples of group analysis for sex

| Canal group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_canal_by_sex.png" width="900"> |

| Discs group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_discs_by_sex.png" width="900"> |

| Foramens group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_foramens_by_sex.png" width="900"> |

| Vertebrae group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_vertebrae_by_sex.png" width="900"> |

### Examples of group analysis for age

| Canal group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_canal_by_age.png" width="900"> |

| Discs group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_discs_by_age.png" width="900"> |

| Foramens group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_foramens_by_age.png" width="900"> |

| Vertebrae group analysis |
| :---: |
| <img src="spinereport/resources/imgs/group_analysis_vertebrae_by_age.png" width="900"> |


