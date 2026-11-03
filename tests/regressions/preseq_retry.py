#!/usr/bin/env python3
"""Check PRESEQ retry argv and recovery using the real module and a stand-in."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--nextflow', default='nextflow')
    parser.add_argument('--nextflow-jar', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = (args.repo/'modules/nf-core/preseq/lcextrap/main.nf').read_text()
    root = (args.output or Path(tempfile.mkdtemp(prefix='atacseq-preseq-retry-'))).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise ValueError('Use an empty output directory')
    command = (['java', '-jar', str(args.nextflow_jar.resolve())] if args.nextflow_jar
               else [shutil.which(args.nextflow) or args.nextflow])
    bins = root/'bin'
    bins.mkdir()
    executable = bins/'preseq'
    executable.write_text('#!'+sys.executable+'\n'+'''import json, os, sys
from pathlib import Path
argv = sys.argv[1:]
if not argv:
    print('preseq Version: 3.2.0', file=sys.stderr)
    raise SystemExit(0)
with open(os.environ['PRESEQ_ARGV_RECORD'], 'a') as stream:
    stream.write(json.dumps(argv)+'\\n')
if not all(option in argv for option in ['-verbose', '-bam', '-seed', '1']):
    raise SystemExit(43)
if os.environ['PRESEQ_FAIL_FIRST'] == '1' and '-defects' not in argv:
    print('Fixture: simulate an initial complexity-model failure', file=sys.stderr)
    raise SystemExit(42)
Path(argv[argv.index('-output')+1]).write_text('Fixture estimate; no scientific computation\\n')
''')
    executable.chmod(0o755)
    env = {key: os.environ[key] for key in ('HOME', 'PATH', 'JAVA_HOME', 'LANG', 'NXF_VER', 'NXF_HOME') if key in os.environ}
    env.update(PATH=str(bins)+':'+env.get('PATH', '/usr/bin:/bin'), NXF_OFFLINE='true',
               NXF_ANSI_LOG='false', NXF_DISABLE_CHECK_LATEST='true')
    results = {}
    cases = [('normal', '0', '0'), ('retry', '1', '0'),
             ('lost_options_control', '1', '1')]
    for name, fail_first, drop_options in cases:
        case = root/name
        case.mkdir()
        (case/'module.nf').write_text(source)
        (case/'input.bam').write_text('Fixture: only the wrapper argv contract is tested\n')
        (case/'main.nf').write_text("include { PRESEQ_LCEXTRAP } from './module'\nworkflow { PRESEQ_LCEXTRAP(Channel.of([[id:'library', single_end:false],file('input.bam')])) }\n")
        (case/'nextflow.config').write_text("process { executor='local'; cpus=1; memory='1 GB'; maxRetries=1; maxErrors=-1; errorStrategy='finish'; ext.args={ System.getenv('PRESEQ_DROP_RETRY_OPTIONS') == '1' && task.attempt > 1 ? '' : '-verbose -bam -seed 1' }; withLabel:error_ignore { errorStrategy='ignore' }; withLabel:error_retry { errorStrategy='retry'; maxRetries=2 } }\nconda.enabled=false\ndocker.enabled=false\nnextflow.enable.telemetry=false\n")
        record = case/'argv.jsonl'
        proc = subprocess.run(command+['run', 'main.nf', '-ansi-log', 'false', '-with-trace', 'trace.tsv', '-work-dir', 'work'],
                              cwd=case, env={**env, 'PRESEQ_ARGV_RECORD': str(record), 'PRESEQ_FAIL_FIRST': fail_first, 'PRESEQ_DROP_RETRY_OPTIONS': drop_options},
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        (case/'run.out').write_bytes(proc.stdout)
        invocations = [json.loads(line) for line in record.read_text().splitlines()] if record.exists() else []
        results[name] = {'exit_code': proc.returncode, 'argv': invocations}
    normal, retry, control = [results[key] for key in ('normal', 'retry', 'lost_options_control')]
    gates = {
        'normal_succeeds_once_without_defects': normal['exit_code'] == 0 and len(normal['argv']) == 1 and '-defects' not in normal['argv'][0],
        'retry_succeeds_with_original_options': retry['exit_code'] == 0 and len(retry['argv']) == 2 and all('-bam' in argv and '-seed' in argv for argv in retry['argv']),
        'defects_only_added_on_retry': len(retry['argv']) == 2 and '-defects' not in retry['argv'][0] and '-defects' in retry['argv'][1],
        'negative_control_detects_lost_options': control['exit_code'] != 0 and len(control['argv']) == 3 and all('-bam' not in argv for argv in control['argv'][1:]),
    }
    result = {'gates': gates, 'cases': results, 'module_sha256': hashlib.sha256(source.encode()).hexdigest(),
              'scope': 'Real Nextflow module; fake PRESEQ executable; no estimate or full pipeline is validated'}
    (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'passed': all(gates.values()), 'gates': gates, 'evidence': str(root/'result.json')}))
    return 0 if all(gates.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
