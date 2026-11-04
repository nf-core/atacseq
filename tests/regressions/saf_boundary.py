#!/usr/bin/env python3
"""Run the actual SAF converter and count a synthetic coordinate-boundary witness."""
import argparse
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
SOURCE = 'modules/local/macs3_consensus.nf'
SUBREAD = 'quay.io/biocontainers/subread@sha256:ccee1f6ebb924fd0b3b6db646a51dabc905697aa25be132386dd490c2318286c'
SAMTOOLS = 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
INTERVALS = [('chrPoc', 0, 10, 'Interval_zero'), ('chrPoc', 100, 110, 'Interval_nonzero')]
BOOLEAN = 'chr\tstart\tend\tpeak_id\tpresent\n' + ''.join(f'{chrom}\t{start}\t{end}\t{name}\t1\n' for chrom, start, end, name in INTERVALS)
BED = ''.join(f'{chrom}\t{start}\t{end}\t{name}\t0\t+\n' for chrom, start, end, name in INTERVALS)
SAM = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrPoc\tLN:1000\n' + ''.join(
    f'{name}\t0\tchrPoc\t{position}\t60\t1M\t*\t0\t0\tA\tI\n'
    for name, position in [('zero_inside', 1), ('edge', 100), ('positive_inside', 105)])


class InfrastructureError(RuntimeError):
    pass


class ScientificFailure(ValueError):
    pass


def digest(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def extract_converter(text):
    lines = text.splitlines()
    headers = [i for i, line in enumerate(lines) if line.strip().startswith('echo -e "GeneID') and line.strip().endswith('> ${prefix}.saf')]
    if len(headers) != 1 or headers[0] + 1 >= len(lines):
        raise InfrastructureError('Unsupported or ambiguous actual SAF header')
    first = headers[0]; awk = lines[first + 1].strip()
    if not (awk.startswith("awk -v FS=") and '${prefix}.boolean.txt >> ${prefix}.saf' in awk and 'FNR > 1' in awk):
        raise InfrastructureError('Unsupported actual SAF conversion command')
    exact = lines[first].strip() + '\n' + awk + '\n'
    # Render only the Nextflow GString escapes and its fixed output-prefix binding.
    rendered = exact.replace('\\$', '$').replace('${prefix}', 'fixture')
    if '${' in rendered:
        raise InfrastructureError('Unresolved production interpolation')
    return exact, rendered


def check_saf(text):
    try:
        rows = list(csv.reader(io.StringIO(text), delimiter='\t', strict=True))
    except csv.Error as error:
        raise ScientificFailure('Malformed actual SAF output') from error
    expected = [['GeneID', 'Chr', 'Start', 'End', 'Strand']] + [
        [name, chrom, str(start + 1), str(end), '+'] for chrom, start, end, name in INTERVALS]
    checks = {'exact_header': bool(rows) and rows[0] == expected[0],
              'complete_unique_ordered_feature_ids': len(rows) == 3 and [r[0] if r else None for r in rows[1:]] == [r[0] for r in expected[1:]],
              'all_start_coordinates_are_BED_start_plus_one': len(rows) == 3 and all(len(r) == 5 and r[2] == e[2] for r, e in zip(rows[1:], expected[1:])),
              'chromosome_end_strand_preserved': len(rows) == 3 and all(len(r) == 5 and [r[1], r[3], r[4]] == [e[1], e[3], e[4]] for r, e in zip(rows[1:], expected[1:]))}
    return {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks, 'expected': expected, 'observed': rows}


def read_counts(text):
    try:
        rows = list(csv.reader([line for line in text.splitlines() if line and not line.startswith('#')], delimiter='\t', strict=True))
    except csv.Error as error:
        raise ScientificFailure('Malformed actual featureCounts output') from error
    if not rows or rows[0][:6] != ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] or len(rows[0]) != 7:
        raise ScientificFailure('Unexpected actual featureCounts output schema')
    result = {}
    for row in rows[1:]:
        if len(row) != 7 or row[0] in result or not re.fullmatch('[0-9]+', row[6]):
            raise ScientificFailure('Malformed or duplicate actual feature count')
        result[row[0]] = int(row[6])
    return result


def produced_text(path):
    try:
        return path.read_text()
    except (FileNotFoundError, UnicodeError) as error:
        raise ScientificFailure('Missing or non-text output after successful tool execution: ' + path.name) from error


def exact_bam_records(text):
    expected = sorted(line for line in SAM.splitlines() if line and not line.startswith('@'))
    observed = sorted(line for line in text.splitlines() if line and not line.startswith('@'))
    return {'passed': observed == expected, 'expected_records': len(expected), 'observed_records': len(observed),
            'all_three_named_alignments_and_fields_preserved': observed == expected}


def count_checks(counts):
    expected = {'Interval_zero': 1, 'Interval_nonzero': 1}
    return {'status': 'PASS' if counts == expected else 'FAIL', 'expected': expected, 'observed': counts,
            'witness': 'The 1M edge read at SAM POS100 belongs to legacy [100,110], but not corrected [101,110]. Positive reads at POS1 and105 remain assigned.'}


def run(args):
    output = args.output
    if output.is_symlink() or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise ValueError('Output must be a new or empty regular directory')
    output.mkdir(parents=True, exist_ok=True); output = output.resolve()
    result = {'schema': 'atacseq-saf-boundary-regression/v1', 'status': 'ERROR', 'commands': [],
              'scope': 'Two constructed intervals and three constructed single-end 1M alignments; actual extracted production SAF command and real Subread counting. No alignment, peak calling, cohort or full pipeline execution.',
              'snapshots_read_or_updated': False, 'probe_source': digest(Path(__file__).read_bytes())}
    def command(argv, label, allowed=(0,)):
        try:
            proc = subprocess.run(argv, cwd=output, capture_output=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired) as error:
            raise InfrastructureError(label + ': ' + type(error).__name__) from error
        for stream in ('stdout', 'stderr'):
            (output / (label + '.' + stream)).write_bytes(getattr(proc, stream))
        result['commands'].append({'label': label, 'argv': list(map(str, argv)), 'exit': proc.returncode,
                                   'stdout': digest(proc.stdout), 'stderr': digest(proc.stderr)})
        if proc.returncode not in allowed:
            raise InfrastructureError(label + ' failed; retained logs identify tool/setup errors')
        return proc.stdout.decode(), proc.stderr.decode()
    try:
        repo = args.repo.resolve(); commit = command(['git', '-C', str(repo), 'rev-parse', args.source_ref or 'HEAD'], 'git-identity')[0].strip()
        tracked = subprocess.check_output(['git', '-C', str(repo), 'show', commit + ':' + SOURCE])
        source = tracked if args.source_ref else (repo / SOURCE).read_bytes()
        (output / 'production-macs3_consensus.nf').write_bytes(source)
        exact, rendered = extract_converter(source.decode())
        (output / 'actual-converter.fragment').write_text(exact)
        (output / 'convert.sh').write_text('set -euo pipefail\n' + rendered)
        result['source_commit'] = commit
        result['source_files'] = [{'path': SOURCE, **digest(source), 'matches_recorded_git_commit': source == tracked}]
        result['extracted'] = {'actual_fragment': digest(exact.encode()), 'rendered_command': digest(rendered.encode()),
                               'only_rendering': 'Nextflow escaped dollar signs and prefix=fixture; awk fields/expression and header unchanged.'}
        for name, data in [('fixture.boolean.txt', BOOLEAN), ('fixture.bed', BED), ('tiny.sam', SAM)]:
            (output / name).write_text(data)
        result['inputs'] = [{'path': name, **digest((output / name).read_bytes())} for name in ['fixture.boolean.txt', 'fixture.bed', 'tiny.sam']]
        result['containers'] = {}
        images = {'subread': args.subread_image, 'samtools': args.samtools_image}
        for name, image in images.items():
            metadata = command(['docker', '--context', args.docker_context, 'image', 'inspect', '--format',
                '{"image_id":{{json .Id}},"repo_digests":{{json .RepoDigests}},"architecture":{{json .Architecture}},"os":{{json .Os}}}', image], 'image-' + name)[0]
            actual = json.loads(metadata)
            if not actual.get('repo_digests') or actual.get('architecture') != 'amd64' or actual.get('os') != 'linux':
                raise InfrastructureError('Missing immutable Linux/amd64 image identity')
            if '@sha256:' in image and not any(x.endswith(image.split('@', 1)[1]) for x in actual['repo_digests']):
                raise InfrastructureError('Requested image digest differs from inspected image')
            result['containers'][name] = {'requested': image, **actual}
        def tool(name, argv, label, entrypoint=None):
            return command(['docker', '--context', args.docker_context, 'run', '--pull', 'never', '--network', 'none',
                '--platform', 'linux/amd64', '--cpus', '1', '--memory', '512m', '--volume', str(output) + ':/probe',
                '--workdir', '/probe', '--entrypoint', entrypoint or argv[0], images[name], *argv[1:]], label)
        stdout, stderr = tool('subread', ['featureCounts', '-v'], 'subread-version')
        if not re.search(r'featureCounts v2\.0\.1\b', stdout + stderr):
            raise InfrastructureError('Expected actual Subread featureCounts2.0.1')
        samtools_version = tool('samtools', ['samtools', '--version'], 'samtools-version')[0].splitlines()[0]
        if not re.fullmatch(r'samtools 1\.24(?:\..*)?', samtools_version):
            raise InfrastructureError('Expected actual samtools1.24')
        result['runtime'] = {'featureCounts': '2.0.1', 'samtools': samtools_version, 'docker_context': args.docker_context}
        tool('subread', ['bash', 'convert.sh'], 'actual-saf-converter')
        saf = output / 'fixture.saf'; result['annotation'] = check_saf(produced_text(saf))
        tool('samtools', ['samtools', 'view', '-b', '-o', 'tiny.bam', 'tiny.sam'], 'construct-bam')
        decoded = tool('samtools', ['samtools', 'view', 'tiny.bam'], 'verify-bam-membership')[0]
        result['counting_input_membership'] = exact_bam_records(decoded)
        if not result['counting_input_membership']['passed']:
            raise InfrastructureError('BAM construction failed exact alignment membership/field check')
        tool('subread', ['featureCounts', '-F', 'SAF', '-a', 'fixture.saf', '-o', 'counts.tsv', '-T', '1', 'tiny.bam'], 'count-boundary')
        result['counting'] = count_checks(read_counts(produced_text(output / 'counts.tsv')))
        result['outputs'] = [{'path': p.name, **digest(p.read_bytes())} for p in [saf, output / 'tiny.bam', output / 'counts.tsv', output / 'counts.tsv.summary']]
        result['status'] = 'PASS' if result['annotation']['status'] == result['counting']['status'] == 'PASS' else 'FAIL'
        result['input_and_source_bytes_unchanged'] = all(digest((output / row['path']).read_bytes())['sha256'] == row['sha256'] for row in result['inputs']) and digest((output / 'production-macs3_consensus.nf').read_bytes()) == digest(source)
        if not result['input_and_source_bytes_unchanged']: raise InfrastructureError('Source/input drift')
    except ScientificFailure as error:
        result.update(status='FAIL', failure_class='scientific_contract', reason=str(error))
    except Exception as error:
        result.update(status='ERROR', reason=type(error).__name__ + ': ' + str(error))
    result['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (output / 'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2]); parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--source-ref'); parser.add_argument('--docker-context', default='default')
    parser.add_argument('--subread-image', default=SUBREAD); parser.add_argument('--samtools-image', default=SAMTOOLS)
    args = parser.parse_args(); result = run(args)
    print(json.dumps({'status': result['status'], 'result': str(args.output / 'RESULT.json')}))
    return {'PASS': 0, 'FAIL': 1, 'ERROR': 2}[result['status']]


if __name__ == '__main__': sys.exit(main())
