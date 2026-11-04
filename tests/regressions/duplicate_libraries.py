#!/usr/bin/env python3
"""Check technical-resequencing duplicate membership through the production Picard modules."""
import argparse
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PICARD = 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'
SAMTOOLS = 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'
SOURCE_FILES = ('nextflow.config', 'conf/modules.config', 'modules/nf-core/picard/addorreplacereadgroups/main.nf',
    'modules/nf-core/picard/markduplicates/main.nf', 'subworkflows/nf-core/fastq_align_chromap/main.nf',
    'subworkflows/nf-core/bam_markduplicates_picard/main.nf', 'workflows/atacseq.nf')
RUNS = ('LIB_A_REP1_T1', 'LIB_A_REP1_T2', 'LIB_B_REP1_T1')
NAMES = {
    'a_primary': 'A_KEEP:1:FC:1:1101:100:100',
    'a_duplicate': 'A_DUP:1:FC:1:1101:101:100',
    'a_unique1': 'A_UNIQUE1:1:FC:1:1101:300:100',
    'a_unique2': 'A_UNIQUE2:1:FC:1:1101:500:100',
    'b_independent': 'B_KEEP:1:FC:1:1101:200:200',
}
READ_RUN = {NAMES[name]: RUNS[run] for name, run in
            [('a_primary', 0), ('a_unique1', 0), ('a_duplicate', 1), ('a_unique2', 1), ('b_independent', 2)]}
READ_START = {NAMES[name]: position for name, position in
              [('a_primary', 100), ('a_unique1', 300), ('a_duplicate', 100), ('a_unique2', 500), ('b_independent', 100)]}


class InfrastructureError(RuntimeError):
    pass


class ScientificFailure(RuntimeError):
    pass


def digest(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def pairs(name, start, quality):
    name = NAMES[name]
    return [
        '\t'.join(map(str, [name, 99, 'chrPoc', start, 60, '20M', '=', start + 50, 70, 'A' * 20, quality * 20])),
        '\t'.join(map(str, [name, 147, 'chrPoc', start + 50, 60, '20M', '=', start, -70, 'T' * 20, quality * 20])),
    ]


def make_inputs(folder):
    layout = {RUNS[0]: pairs('a_primary', 100, 'I') + pairs('a_unique1', 300, 'I'),
              RUNS[1]: pairs('a_duplicate', 100, '!') + pairs('a_unique2', 500, 'I'),
              RUNS[2]: pairs('b_independent', 100, '5')}
    for run, rows in layout.items():
        rows.sort(key=lambda row: (int(row.split('\t')[3]), row.split('\t')[0]))
        (folder / (run + '.sam')).write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrPoc\tLN:2000\n' + '\n'.join(rows) + '\n')


def parse_sam(text):
    records, groups = [], {}
    for line in text.splitlines():
        if line.startswith('@RG\t'):
            try:
                fields = dict(field.split(':', 1) for field in line.split('\t')[1:])
                identifier = fields['ID']
            except (ValueError, KeyError) as error:
                raise ScientificFailure('Malformed produced read-group header') from error
            if identifier in groups:
                raise ScientificFailure('Duplicate read-group ID in observed BAM header')
            groups[identifier] = fields
        elif line and not line.startswith('@'):
            fields = line.split('\t')
            if len(fields) < 11:
                raise ScientificFailure('Malformed produced SAM alignment')
            try:
                flag, position = int(fields[1]), int(fields[3])
            except ValueError as error:
                raise ScientificFailure('Malformed produced SAM flag or coordinate') from error
            tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in fields[11:] if x.count(':') >= 2}
            records.append({'name': fields[0], 'mate': 1 if flag & 64 else 2 if flag & 128 else None,
                'flag': flag, 'duplicate': bool(flag & 1024), 'rg': tags.get('RG'),
                'reference': fields[2], 'position': position, 'cigar': fields[5]})
    return groups, records


def scientific_checks(groups, records):
    expected_groups = {run: {'--RGLB': run.rsplit('_T', 1)[0], '--RGSM': run.rsplit('_T', 1)[0]} for run in RUNS}
    expected_names = set(NAMES.values())
    expected_duplicate = {NAMES['a_duplicate']}
    by_name = {name: [r for r in records if r['name'] == name] for name in expected_names}
    observed_names = {r['name'] for r in records}
    observed_duplicate = {r['name'] for r in records if r['duplicate']}
    flags_consistent = all({r['mate'] for r in rows} == {1, 2} and len(rows) == 2
                           and len({r['duplicate'] for r in rows}) == 1 for rows in by_name.values())
    comparisons = {
        'all_input_pairs_present_exactly_once': observed_names == expected_names and len(records) == 10 and flags_consistent,
        'only_lower_quality_technical_resequence_pair_marked_duplicate': observed_duplicate == expected_duplicate,
        'independent_biological_library_same_coordinates_retained': all(not r['duplicate'] for r in by_name[NAMES['b_independent']]) and len(by_name[NAMES['b_independent']]) == 2,
        'retained_pair_membership_exact': {r['name'] for r in records if not r['duplicate']} == expected_names - expected_duplicate,
        'read_group_ids_exact': set(groups) == set(expected_groups),
        'read_group_tags_match_actual_selector_metadata': all(groups.get(k, {}).get('LB') == v['--RGLB'] and groups.get(k, {}).get('SM') == v['--RGSM'] for k, v in expected_groups.items()),
        'every_alignment_has_a_declared_read_group': all(r['rg'] in groups for r in records),
        'every_alignment_keeps_its_actual_run_read_group': all(r['rg'] == READ_RUN.get(r['name']) for r in records),
        'alignment_coordinates_cigar_and_nonduplicate_flags_unchanged': all(
            r['reference'] == 'chrPoc' and r['cigar'] == '20M'
            and r['position'] == READ_START.get(r['name'], -999) + (50 if r['mate'] == 2 else 0)
            and r['flag'] & ~1024 == (99 if r['mate'] == 1 else 147) for r in records),
        'technical_runs_share_library_and_keep_distinct_ids': len(set(expected_groups)) == 3 and groups.get(RUNS[0], {}).get('LB') == groups.get(RUNS[1], {}).get('LB'),
        'biological_library_is_distinct': groups.get(RUNS[0], {}).get('LB') != groups.get(RUNS[2], {}).get('LB'),
    }
    return {'status': 'PASS' if all(comparisons.values()) else 'FAIL', 'checks': comparisons,
        'expected_duplicate_pairs': sorted(expected_duplicate), 'observed_duplicate_pairs': sorted(observed_duplicate),
        'expected_retained_pairs': sorted(expected_names - expected_duplicate),
        'observed_retained_pairs': sorted({r['name'] for r in records if not r['duplicate']}),
        'duplicate_alignment_records': sum(r['duplicate'] for r in records), 'alignment_records': len(records),
        'observed_read_groups': groups, 'observed_alignment_membership': records}


def run(args):
    output = args.output.resolve()
    if output.is_symlink() or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise ValueError('Output must be a new or empty directory')
    output.mkdir(parents=True, exist_ok=True)
    result = {'schema': 'atacseq-duplicate-library-regression/v1', 'status': 'ERROR', 'commands': [],
              'snapshots_read_or_updated': False, 'actual_nextflow_module_wrappers': True}
    docker = ['docker'] + (['--context', args.docker_context] if args.docker_context else [])

    def command(argv, name, env=None, allowed=(0,)):
        try:
            proc = subprocess.run(argv, cwd=output, env=env, capture_output=True, timeout=600, check=False)
        except (OSError, subprocess.TimeoutExpired) as error:
            raise InfrastructureError(name + ': ' + type(error).__name__) from error
        for suffix, data in [('stdout', proc.stdout), ('stderr', proc.stderr)]:
            (output / (name + '.' + suffix)).write_bytes(data)
        result['commands'].append({'name': name, 'argv': list(map(str, argv)), 'exit': proc.returncode,
                                  'stdout': digest(output / (name + '.stdout')), 'stderr': digest(output / (name + '.stderr'))})
        if proc.returncode not in allowed:
            raise InfrastructureError(name + ' failed; see retained stdout/stderr')
        return proc.stdout.decode()

    try:
        repo = args.repo.resolve()
        result['source_commit'] = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
        source = output / 'source'
        result['source_files'] = []
        for relative in SOURCE_FILES:
            target = source / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(repo / relative, target)
            tracked = subprocess.check_output(['git', '-C', str(repo), 'show', result['source_commit'] + ':' + relative])
            result['source_files'].append({'path': relative, **digest(target), 'matches_recorded_git_commit': target.read_bytes() == tracked})
        for name in ('duplicate_libraries.py', 'duplicate_libraries.nf'):
            target = source / 'tests/regressions' / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(HERE / name, target)
        result['guard_files'] = [{'path': 'tests/regressions/' + name, **digest(source / 'tests/regressions' / name)}
                                 for name in ('duplicate_libraries.py', 'duplicate_libraries.nf')]
        for module in ('addorreplacereadgroups', 'markduplicates'):
            text = (source / ('modules/nf-core/picard/' + module + '/main.nf')).read_text()
            declarations = {x for x in re.findall(r"'([^']*picard:[^']+)'", text) if not x.startswith('https://')}
            if declarations != {PICARD}:
                raise InfrastructureError('Unexpected Picard declaration; review the image and fixture explicitly')
        result['containers'] = {}
        for name, image in [('picard', PICARD), ('samtools', SAMTOOLS)]:
            inspected = json.loads(command(docker + ['image', 'inspect', image], name + '-image'))[0]
            info = {key: inspected[field] for key, field in
                    [('id', 'Id'), ('digests', 'RepoDigests'), ('architecture', 'Architecture')]}
            if not info['digests']:
                raise InfrastructureError('Cached image lacks an immutable repository digest')
            result['containers'][name] = dict(info, declaration=image)

        def samtools(argv, name):
            image = result['containers']['samtools']
            return command(docker + ['run', '--rm', '--pull', 'never', '--network', 'none',
                '--platform', 'linux/' + image['architecture'], '--cpus', '1', '--memory', '1g',
                '--mount', 'type=bind,source=' + str(output) + ',target=/evidence', '-w', '/evidence',
                '--entrypoint', 'samtools', image['digests'][0]] + argv, name)

        make_inputs(output)
        result['inputs'] = []
        for run_id in RUNS:
            result['inputs'].append({'path': run_id + '.sam', **digest(output / (run_id + '.sam'))})
            samtools(['view', '-b', '-o', run_id + '.bam', run_id + '.sam'], 'prepare-' + run_id)
        defaults = re.findall(r'(?ms)^params\s*\{.*?^\}', (source / 'nextflow.config').read_text())
        if len(defaults) != 1:
            raise InfrastructureError('Cannot isolate all production parameter defaults')
        config = output / 'guard.config'
        config.write_text(defaults[0] + '\n' + '\n'.join([
            'params.outdir = ' + json.dumps(str(output / 'published')),
            'params.aligner = "chromap"',
            'includeConfig ' + json.dumps(str(source / 'conf/modules.config')),
            'process.executor = "local"', 'process.cpus = 1', 'process.memory = "1 GB"',
            'process.maxForks = 1', 'docker.enabled = true', 'docker.pull = false',
            'docker.runOptions = "--pull never --network none --hostname duplicate-guard --add-host duplicate-guard:127.0.0.1"',
            'singularity.enabled = false', 'conda.enabled = false',
            'trace.fields = "task_id,process,tag,status,exit,hash,workdir,container"',
            'params.samtools_image = ' + json.dumps(SAMTOOLS),
            'params.fixture_dir = ' + json.dumps(str(output)),
            'params.result_path = ' + json.dumps(str(output / 'produced.json'))]) + '\n')
        result['effective_config'] = digest(config)
        env = dict(os.environ)
        if args.docker_context:
            env['DOCKER_CONTEXT'] = args.docker_context
        argv = [args.nextflow, '-C', str(config), 'run', str(source / 'tests/regressions/duplicate_libraries.nf'),
                '-work-dir', str(output / 'work'), '-with-trace', str(output / 'trace.tsv'), '-ansi-log', 'false']
        command(argv, 'nextflow', env=env, allowed=(0, 1))
        trace = output / 'trace.tsv'
        tasks = list(csv.DictReader(trace.open(), delimiter='\t')) if trace.is_file() else []
        result['tasks'] = []
        for task in tasks:
            work = Path(task['workdir'])
            scripts = {name: digest(work / name) for name in ('.command.sh', '.command.run', '.exitcode') if (work / name).is_file()}
            result['tasks'].append(dict(task, files=scripts))
        wrappers = [t for t in tasks if t['process'].endswith((':PICARD_ADDORREPLACEREADGROUPS', ':PICARD_MARKDUPLICATES'))]
        for task in wrappers:
            if task['status'] == 'COMPLETED' and task['exit'] == '0' and not any(
                    p.is_file() and not p.is_symlink() and p.stat().st_size for p in Path(task['workdir']).glob('*.bam')):
                raise ScientificFailure('Successful Picard module did not produce a nonempty BAM')
        qualified = 'NFCORE_ATACSEQ:ATACSEQ:'
        if (len(wrappers) != 4 or any(t['status'] != 'COMPLETED' or t['exit'] != '0' for t in wrappers)
                or sum(t['process'] == qualified + 'FASTQ_ALIGN_CHROMAP:PICARD_ADDORREPLACEREADGROUPS' for t in wrappers) != 3
                or sum(t['process'] == qualified + 'MERGED_LIBRARY_MARKDUPLICATES_PICARD:PICARD_MARKDUPLICATES' for t in wrappers) != 1):
            raise InfrastructureError('Expected four successful production-qualified Picard module tasks')
        if result['commands'][-1]['exit'] != 0:
            raise InfrastructureError('Nextflow did not finish successfully')
        if not (output / 'produced.json').is_file():
            raise ScientificFailure('Successful Picard modules did not emit the marked BAM')
        produced = json.loads((output / 'produced.json').read_text())
        result['nextflow_version'] = produced['nextflow_version']
        bam = Path(produced['bam'])
        if not bam.is_file() or bam.stat().st_size == 0:
            raise ScientificFailure('Successful Picard modules produced a missing or empty marked BAM')
        try:
            sam = samtools(['view', '-h', str(bam.relative_to(output))], 'marked-sam')
        except InfrastructureError as error:
            if result['commands'][-1]['exit'] in (1, 2):
                raise ScientificFailure('Successful Picard modules produced an undecodable BAM') from error
            raise
        groups, records = parse_sam(sam)
        result['scientific'] = scientific_checks(groups, records)
        result['status'] = result['scientific']['status']
        result['outputs'] = [{'path': str(p.relative_to(output)), **digest(p)} for p in sorted((output / 'work').glob('*/*/*'))
                             if p.is_file() and not p.is_symlink() and (p.suffix == '.bam' or p.name.endswith('.metrics.txt'))]
        if any(digest(source / row['path'])['sha256'] != row['sha256'] for row in result['source_files']):
            raise InfrastructureError('Frozen production source changed during execution')
    except ScientificFailure as error:
        result.update(status='FAIL', reason=str(error), error_type=type(error).__name__)
    except Exception as error:
        result.update(status='ERROR', reason=str(error), error_type=type(error).__name__)
    result['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write_json(output / 'RESULT.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=HERE.parents[1])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--docker-context', default='')
    parser.add_argument('--nextflow', default='nextflow')
    args = parser.parse_args()
    try:
        result = run(args)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps({'status': result['status'], 'result': str(args.output / 'RESULT.json')}))
    return {'PASS': 0, 'FAIL': 1, 'ERROR': 2}[result['status']]


if __name__ == '__main__':
    sys.exit(main())
