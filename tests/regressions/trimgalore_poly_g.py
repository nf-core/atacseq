#!/usr/bin/env python3
"""Retain legitimate poly-G/poly-C reads through the actual Trim Galore module."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

SEQUENCES = (
    'TACCTGACTGTCAGTACGATGCTACCTAGTCACTGACCTGATCGTACGACTACGTC' + 'G' * 24,
    'C' * 24 + 'CGTATGACCTACAGTCGATCCTAGTCAGTTCGACATGCTGACCTAGTACGTCAGTC',
)
SOURCES = ('conf/modules.config', 'modules/nf-core/trimgalore/main.nf',
           'modules/nf-core/trimgalore/environment.yml', 'nextflow.config')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def output_identity(path):
    result = {'path': str(path)}
    try:
        result['sha256'] = sha(path)
    except OSError as error:
        result['read_error'] = str(error)
    return result


def expected_records(mate):
    sequence = SEQUENCES[mate - 1]
    return [(f'poly_g_{index:03}/{mate}', sequence, 'I' * len(sequence)) for index in range(32)]


def write_reads(path, mate):
    text = ''.join(f'@{name}\n{sequence}\n+\n{quality}\n' for name, sequence, quality in expected_records(mate))
    Path(path).write_bytes(gzip.compress(text.encode(), mtime=0))


def read_records(path):
    with gzip.open(path, 'rt', encoding='ascii') as stream:
        lines = stream.read().splitlines()
    if len(lines) % 4:
        raise ValueError('Incomplete FASTQ record')
    rows = []
    for start in range(0, len(lines), 4):
        name, sequence, plus, quality = lines[start:start + 4]
        if not name.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
            raise ValueError('Malformed FASTQ record')
        rows.append((name[1:], sequence, quality))
    return rows


def assess(path, mate):
    expected = expected_records(mate)
    try:
        observed = read_records(path)
        return {'passed': Counter(expected) == Counter(observed), 'expected': expected, 'observed': observed}
    except (ValueError, OSError, UnicodeError, EOFError) as error:
        return {'passed': False, 'expected': expected, 'error': str(error)}


def source_identity(repo):
    head = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    dirty = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain'], capture_output=True, text=True)
    return {'files': {name: sha(repo / name) for name in SOURCES},
            'observed_head': head.stdout.strip() if head.returncode == 0 else None,
            'clean_commit': head.stdout.strip() if head.returncode == dirty.returncode == 0 and not dirty.stdout else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--docker-context', default='default')
    parser.add_argument('--nextflow', default='nextflow')
    parser.add_argument('--nextflow-version', default=os.environ.get('NXF_VER', '25.10.4'))
    parser.add_argument('--java-home', type=Path)
    args = parser.parse_args()
    repo, root = args.repo.resolve(), args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        parser.error('--output must be empty')
    result = {'status': 'ERROR', 'commands': [], 'scope': 'Actual module/config/defaults and declared image; 32 synthetic pairs. No snapshots, full pipeline or biological estimate is evaluated.'}

    def run(command, label, env=None):
        record = {'command': command}
        result['commands'].append(record)
        process = subprocess.run(command, cwd=root, env=env, capture_output=True, timeout=240)
        for stream in ('stdout', 'stderr'):
            path = root / f'{label}.{stream}.log'
            path.write_bytes(getattr(process, stream))
            record[stream] = {'file': path.name, 'sha256': sha(path)}
        record['exit'] = process.returncode
        if process.returncode:
            raise ValueError(f'{label} execution failed; inspect retained logs')
        return process.stdout.decode()

    try:
        result['source_before'] = source_identity(repo)
        harness = [Path(__file__), Path(__file__).with_name('trim_poly_g.nf')]
        result['harness_before'] = {path.name: sha(path) for path in harness}
        module = (repo / SOURCES[1]).read_text()
        images = set(value for value in re.findall(r"'([^']*trim-galore:[^']+)'", module) if not value.startswith('https://'))
        if len(images) != 1:
            raise ValueError('Expected one source-declared Docker image')
        image = images.pop()
        docker = ['docker', '--context', args.docker_context]
        inspect = docker + ['image', 'inspect', '--format', '{"id":{{json .Id}},"repo_digests":{{json .RepoDigests}},"architecture":{{json .Architecture}}}', image]
        result['container'] = {'declared_reference': image, **json.loads(run(inspect, 'image-before'))}
        result['tool_identity'] = run(docker + ['run', '--rm', '--platform', 'linux/amd64', '--network', 'none', '--cpus', '1', '--memory', '2g', '--entrypoint', 'sh', image, '-c', 'command -v trim_galore; sha256sum "$(command -v trim_galore)"; trim_galore --version'], 'tool')
        inputs = [root / f'input_{mate}.fastq.gz' for mate in (1, 2)]
        for mate, path in enumerate(inputs, 1):
            write_reads(path, mate)
        result['inputs'] = {path.name: sha(path) for path in inputs}
        defaults = re.findall(r'(?ms)^params\s*\{.*?^\}', (repo / 'nextflow.config').read_text())
        if len(defaults) != 1:
            raise ValueError('Cannot isolate production parameter defaults')
        config = root / 'probe.config'
        config.write_text(defaults[0] + '\nparams.outdir = ' + json.dumps(str(root / 'published')) + '\nincludeConfig ' + json.dumps(str(repo / 'conf/modules.config')) + '''
process.executor = 'local'
process.cpus = 1
process.memory = '2 GB'
process.time = '2 min'
process.maxForks = 1
conda.enabled = false
docker.enabled = true
docker.runOptions = '--platform linux/amd64 --network none --cpus 1 --memory 2g'
nextflow.enable.telemetry = false
trace.fields = 'task_id,hash,name,status,exit,container,workdir,attempt'
''')
        result['config_sha256'] = sha(config)
        env = {**os.environ, 'DOCKER_CONTEXT': args.docker_context, 'NXF_VER': args.nextflow_version,
               'NXF_HOME': str(root / 'nxf-home'), 'NXF_DIST': os.environ.get('NXF_DIST', str(Path.home() / '.nextflow/framework')),
               'NXF_OFFLINE': 'true', 'NXF_ANSI_LOG': 'false', 'NXF_DISABLE_CHECK_LATEST': 'true', 'NXF_JVM_ARGS': '-Xms128m -Xmx1g'}
        if args.java_home:
            env['JAVA_HOME'] = str(args.java_home.resolve())
        result['nextflow_version'] = run([args.nextflow, '-version'], 'nextflow-version', env)
        if f'version {args.nextflow_version}' not in result['nextflow_version']:
            raise ValueError('Actual Nextflow version differs from the requested version')
        receipt = root / 'reads.json'
        run([args.nextflow, '-C', str(config), 'run', str(harness[1]), '-ansi-log', 'false', '-with-trace', 'trace.tsv', '-work-dir', str(root / 'work'), '--repo', str(repo), '--read1', str(inputs[0]), '--read2', str(inputs[1]), '--receipt', str(receipt)], 'nextflow-run', env)
        with (root / 'trace.tsv').open() as stream:
            tasks = list(csv.DictReader(stream, delimiter='\t'))
        if len(tasks) != 1 or tasks[0]['status'] != 'COMPLETED' or tasks[0]['exit'] != '0' or tasks[0]['container'] != image:
            raise ValueError('Expected one completed uncached task using the declared image')
        result['task'] = tasks[0]
        work = Path(tasks[0]['workdir'])
        result['actual_command'] = (work / '.command.sh').read_text()
        result['task_file_sha256'] = {name: sha(work / name) for name in ('.command.sh', '.command.run', '.command.err', '.command.out', '.exitcode')}
        outputs = [Path(path) for path in json.loads(receipt.read_text())['reads']]
        result['mates'] = []
        for mate in (1, 2):
            paths = [path for path in outputs if path.name.startswith(f'poly_g_{mate}_')]
            outcome = assess(paths[0], mate) if len(paths) == 1 else {'passed': False, 'reason': 'missing_or_ambiguous_mate', 'expected': expected_records(mate)}
            result['mates'].append({'mate': mate, 'outputs': [output_identity(path) for path in paths], **outcome})
        result['source_after'] = source_identity(repo)
        result['harness_after'] = {path.name: sha(path) for path in harness}
        after = json.loads(run(inspect, 'image-after'))
        if result['source_before']['files'] != result['source_after']['files'] or result['harness_before'] != result['harness_after'] or after['id'] != result['container']['id']:
            raise ValueError('Source, harness or image changed during execution')
        result['status'] = 'PASS' if all(mate['passed'] for mate in result['mates']) else 'FAIL'
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        result['status'], result['error'] = 'ERROR', f'{type(error).__name__}: {error}'
    (root / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'result': str(root / 'RESULT.json')}))
    return {'PASS': 0, 'FAIL': 1, 'ERROR': 2}[result['status']]


if __name__ == '__main__':
    raise SystemExit(main())
