#!/usr/bin/env python3
"""Check pins/patches; run in a Python environment containing nf-core/tools."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import shutil
from urllib.request import urlopen

import nf_core
from nf_core.components.components_differ import ComponentsDiffer


MODULE_FILES = {
    'gffread': ['main.nf', 'environment.yml', 'tests/main.nf.test.snap'],
    'preseq/lcextrap': ['main.nf'],
    'trimgalore': ['main.nf', 'environment.yml', 'tests/main.nf.test.snap'],
}


def regenerated_replay(name, downloaded, patched, output):
    """Exercise the patch regeneration performed during nf-core updates."""
    installed = output / 'patched' / name
    shutil.copytree(downloaded, installed, dirs_exist_ok=True)
    for relative, lines in patched.items():
        target = installed / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(''.join(lines))
    patch = output / 'regenerated' / (name.replace('/', '-') + '.diff')
    patch.parent.mkdir(parents=True, exist_ok=True)
    component = Path('modules/nf-core') / name
    ComponentsDiffer.write_diff_file(
        patch, name, 'nf-core', downloaded, installed, file_action='w',
        for_git=False, dsp_from_dir=component, dsp_to_dir=component
    )
    forward = ComponentsDiffer.try_apply_patch('modules', name, 'nf-core', patch, downloaded)
    reverse = ComponentsDiffer.try_apply_patch('modules', name, 'nf-core', patch, installed, reverse=True)
    return {
        'forward': {relative: ''.join(forward[relative]).encode() == (installed / relative).read_bytes()
                    for relative in patched},
        'reverse': {relative: ''.join(reverse[relative]).encode() == (downloaded / relative).read_bytes()
                    for relative in patched},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--cache', type=Path, help='Reuse downloaded upstream files')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = args.output or Path(tempfile.mkdtemp(prefix='atacseq-module-patches-'))
    output.mkdir(parents=True, exist_ok=True)
    cache = args.cache or output / 'upstream'
    registry = json.loads((args.repo / 'modules.json').read_text())
    modules = registry['repos']['https://github.com/nf-core/modules.git']['modules']['nf-core']
    results = {}
    for name, files in MODULE_FILES.items():
        entry = modules[name]
        component = Path('modules/nf-core') / name
        downloaded = cache / entry['git_sha'] / component
        hashes = {}
        for relative in files:
            target = downloaded / relative
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                url = f"https://raw.githubusercontent.com/nf-core/modules/{entry['git_sha']}/{component / relative}"
                with urlopen(url, timeout=30) as response:
                    target.write_bytes(response.read())
            hashes[relative] = hashlib.sha256(target.read_bytes()).hexdigest()
        stage = 'registered patch replay'
        try:
            patched = (ComponentsDiffer.try_apply_patch(
                'modules', name, 'nf-core', args.repo / entry['patch'], downloaded
            ) if entry.get('patch') else {})
            matches = {
                relative: (''.join(patched[relative]).encode() if relative in patched
                           else (downloaded / relative).read_bytes())
                == (args.repo / component / relative).read_bytes()
                for relative in files
            }
            stage = 'canonical regeneration replay'
            regeneration = (regenerated_replay(name, downloaded, patched, output)
                            if patched else {})
            results[name] = {'passed': all(matches.values()) and all(
                                 all(checks.values()) for checks in regeneration.values()),
                             'matches': matches, 'regeneration': regeneration,
                             'git_sha': entry['git_sha'], 'upstream_sha256': hashes}
        except LookupError as error:
            results[name] = {'passed': False, 'stage': stage, 'error': str(error),
                             'git_sha': entry['git_sha'], 'upstream_sha256': hashes}
    record = output / 'result.json'
    record.write_text(json.dumps({'nf_core_version': nf_core.__version__, 'components': results}, indent=2) + '\n')
    passed = all(result['passed'] for result in results.values())
    print(json.dumps({'passed': passed, 'nf_core_version': nf_core.__version__,
                      'components': results, 'evidence': str(record)}))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
