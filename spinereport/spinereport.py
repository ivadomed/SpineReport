import argparse
import multiprocessing as mp
from pathlib import Path
import textwrap
import os, importlib

import totalspineseg.resources as resources

from spinereport.utils.generate_reports import generate_reports
from spinereport.utils.measure_seg import measure_seg_mp


def _add_group_seg_args(parser, group):
    '''Add the per-group segmentation-source flags (--{group}-*-dir).'''
    parser.add_argument(
        f'--{group}-tss', f'-{group[0]}', type=Path, default=None,
        help=f'The output folder where totalspineseg outputs are located for the {group} group.'
    )
    parser.add_argument(
        f'--{group}-images-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group NIfTI images resampled to 1 mm isotropic. If you don"t use --{group}-tss'
    )
    parser.add_argument(
        f'--{group}-labels-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group NIfTI landmark labels (posterior tip of the discs). Use only to replace totalspineseg labels.'
    )
    parser.add_argument(
        f'--{group}-sc-seg-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group BINARY spinal-cord segmentations. Use only to replace totalspineseg SC seg.'
    )
    parser.add_argument(
        f'--{group}-canal-seg-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group BINARY spinal-canal segmentations (SC + CSF). Use only to replace totalspineseg canal seg.'
    )
    parser.add_argument(
        f'--{group}-vertebrae-seg-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group multi-label vertebrae segmentations. Must contain a map.json (see README.md). Use only to replace totalspineseg vertebrae seg.'
    )
    parser.add_argument(
        f'--{group}-discs-seg-dir', type=Path, default=None,
        help=f'Flat folder of {group}-group multi-label discs segmentations. Must contain a map.json (see README.md). Use only to replace totalspineseg discs seg.'
    )


def _resolve_group_paths(group, args):
    '''
    Resolve the effective per-structure paths for a group.

    Individual --{group}-*-dir flags always win. Anything not set falls back to the
    totalspineseg subfolder under --{group}-dir when that flag is provided.
    Returns a dict with keys: labels, combined, sc, canal, vertebrae, discs
    (each is a Path or None).
    '''
    root = getattr(args, f'{group}_dir')

    def _pick(explicit, tss_subfolder):
        if explicit is not None:
            return explicit
        if root is not None and tss_subfolder is not None:
            return root / tss_subfolder
        return None

    paths = {
        'images': _pick(getattr(args, f'{group}_images_dir'), 'input'),
        'labels': _pick(getattr(args, f'{group}_labels_dir'), 'step1_levels'),
        'combined': _pick(None, 'step2_output') if root is not None else None,
        'sc': getattr(args, f'{group}_sc_seg_dir'),
        'canal': getattr(args, f'{group}_canal_seg_dir'),
        'vertebrae': getattr(args, f'{group}_vertebrae_seg_dir'),
        'discs': getattr(args, f'{group}_discs_seg_dir'),
    }

    if paths['images'] is None:
        raise ValueError(
            f'{group} group: no image folder specified. Provide --{group}-dir (totalspineseg output) '
            f'or --{group}-images-dir (flat folder of 1 mm isotropic images).'
        )
    if paths['labels'] is None:
        raise ValueError(
            f'{group} group: no label folder specified. Provide --{group}-dir (totalspineseg output) '
            f'or --{group}-labels-dir.'
        )
    if all(paths[k] is None for k in ('combined', 'sc', 'canal', 'vertebrae', 'discs')):
        raise ValueError(
            f'{group} group: no segmentation source specified. Provide --{group}-dir (combined seg via '
            f'step2_output) or at least one of --{group}-sc-seg-dir / --{group}-canal-seg-dir / '
            f'--{group}-vertebrae-seg-dir / --{group}-discs-seg-dir.'
        )
    return paths


def _pick_metrics_root(group_paths, group, ofolder):
    '''Choose where to cache the per-group metrics_output folder.'''
    root = group_paths['combined']
    if root is not None:
        return root.parent / 'metrics_output'
    # No totalspineseg folder: stash cached metrics under the report ofolder.
    return Path(ofolder) / f'{group}_metrics_output'


def main():
    # Description and arguments
    parser = argparse.ArgumentParser(
        description=' '.join(f'''
            This script processes segmentation folders to generate reports.
            Segmentations can come from a totalspineseg output folder (--test-dir / --control-dir) or from
            separate per-structure folders (--*-sc-seg-dir, --*-canal-seg-dir, --*-vertebrae-seg-dir,
            --*-discs-seg-dir); the two can also be mixed, with per-structure folders overriding what is
            in the totalspineseg folder. totalspineseg can be skipped entirely when every structure is
            provided along with --*-images-dir and --*-labels-dir, as long as all inputs share the same
            1 mm isotropic space. All folders must be flat (no per-subject subdirectories) and filenames
            must follow the BIDS naming convention.
        '''.split()),
        epilog=textwrap.dedent('''
            Examples:
            # totalspineseg only
            spinereport -t test_totalspineseg -c control_totalspineseg -o reports

            # totalspineseg with custom SC segmentation
            spinereport \\
                --test-sc-seg-dir test/sc --control-sc-seg-dir ctrl/sc \\
                -t test_totalspineseg -c control_totalspineseg -o reports

        '''),
        formatter_class=argparse.RawTextHelpFormatter
    )
    _add_group_seg_args(parser.add_argument_group('test group'), 'test')
    _add_group_seg_args(parser.add_argument_group('control group'), 'control')
    parser.add_argument(
        '--ofolder', '-o', type=Path, required=True,
        help='The folder where reports will be saved (required).'
    )
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Whether to overwrite existing folders.'
    )
    parser.add_argument(
        '--max-workers', '-w', type=int, default=mp.cpu_count(),
        help='Max worker to run in parallel proccess, defaults to multiprocessing.cpu_count().'
    )
    parser.add_argument(
        '--quiet', '-q', action="store_true", default=False,
        help='Do not display inputs and progress bar, defaults to false (display).'
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    ofolder = args.ofolder
    overwrite = args.overwrite
    max_workers = args.max_workers
    quiet = args.quiet

    # Resolve effective per-group paths (per-structure flags win over --{group}-dir defaults).
    test_paths = _resolve_group_paths('test', args)
    control_paths = _resolve_group_paths('control', args)

    # Use default mapping path
    resources_path = importlib.resources.files(resources)
    mapping_path = os.path.join(resources_path, 'labels_maps/tss_map.json')

    # Print the argument values if not quiet
    if not quiet:
        print(textwrap.dedent(f'''
            Running {Path(__file__).stem} with the following params:
            test_paths = {test_paths}
            control_paths = {control_paths}
            ofolder = "{ofolder}"
            overwrite = "{overwrite}"
            mapping_path = "{mapping_path}"
            max_workers = {max_workers}
            quiet = {quiet}
        '''))

    # Run spinereport
    run_spinereport(
        test_paths=test_paths,
        control_paths=control_paths,
        ofolder=ofolder,
        overwrite=overwrite,
        mapping_path=mapping_path,
        max_workers=max_workers,
        quiet=quiet
    )


def _measure_group(group, paths, ofolder, overwrite, mapping_path, max_workers, quiet):
    '''Run measure_seg_mp for a group unless a cached metrics_output already exists.'''
    for key in ('images', 'labels'):
        if not paths[key].exists():
            raise FileNotFoundError(f'{group} group {key} folder "{paths[key]}" does not exist.')

    metrics_path = _pick_metrics_root(paths, group, ofolder)
    if metrics_path.exists() and not overwrite:
        return metrics_path

    if not quiet: print(f'\nMeasuring segmentations for {group} group in "{paths["images"]}"...')
    measure_seg_mp(
        images_path=paths['images'],
        labels_path=paths['labels'],
        ofolder_path=metrics_path,
        segs_path=paths['combined'],
        sc_segs_path=paths['sc'],
        canal_segs_path=paths['canal'],
        vertebrae_segs_path=paths['vertebrae'],
        discs_segs_path=paths['discs'],
        mapping_path=mapping_path,
        max_workers=max_workers,
        quiet=quiet,
    )
    return metrics_path


def run_spinereport(
        test_paths: dict,
        control_paths: dict,
        ofolder: Path,
        overwrite: bool,
        mapping_path: str,
        max_workers: int,
        quiet: bool
    ):

    test_metrics_path = _measure_group('test', test_paths, ofolder, overwrite, mapping_path, max_workers, quiet)

    # If the two groups share the same measurement inputs (e.g. same totalspineseg folder), reuse.
    if control_paths == test_paths:
        control_metrics_path = test_metrics_path
    else:
        control_metrics_path = _measure_group('control', control_paths, ofolder, overwrite, mapping_path,
                                              max_workers, quiet)

    # Generate reports
    if not quiet: print(f'\nGenerating reports in "{ofolder}"...')
    generate_reports(
        test_path=test_metrics_path,
        control_path=control_metrics_path,
        ofolder_path=ofolder,
        max_workers=max_workers,
        quiet=quiet
    )

    if not quiet: print('Reports generation completed.')
