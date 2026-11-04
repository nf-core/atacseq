import importlib.util
from pathlib import Path
import unittest
import tempfile

spec = importlib.util.spec_from_file_location('saf_probe', Path(__file__).with_name('saf_boundary.py'))
probe = importlib.util.module_from_spec(spec); spec.loader.exec_module(probe)


class SafBoundaryContracts(unittest.TestCase):
    def setUp(self):
        self.valid = 'GeneID\tChr\tStart\tEnd\tStrand\nInterval_zero\tchrPoc\t1\t10\t+\nInterval_nonzero\tchrPoc\t101\t110\t+\n'

    def test_zero_and_nonzero_start_are_independently_one_based(self):
        self.assertEqual(probe.check_saf(self.valid)['status'], 'PASS')
        for text in [self.valid.replace('\t1\t10\t', '\t0\t10\t'), self.valid.replace('\t101\t110\t', '\t100\t110\t')]:
            self.assertEqual(probe.check_saf(text)['status'], 'FAIL')

    def test_dropped_duplicate_swapped_or_changed_annotation_rejects(self):
        rows = self.valid.splitlines(True)
        for text in [''.join(rows[:-1]), ''.join(rows + rows[-1:]), ''.join([rows[0], rows[2], rows[1]]),
                     self.valid.replace('chrPoc', 'chrOther'), self.valid.replace('\t110\t', '\t111\t'), self.valid.replace('\t+\n', '\t-\n')]:
            self.assertEqual(probe.check_saf(text)['status'], 'FAIL')

    def test_actual_converter_is_rendered_without_rewriting_awk(self):
        actual = '''    echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > ${prefix}.saf
    awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 { print \\$4, \\$1, \\$2+1, \\$3, "+" }' ${prefix}.boolean.txt >> ${prefix}.saf
'''
        exact, rendered = probe.extract_converter(actual)
        self.assertIn('\\$2+1', exact); self.assertIn('$2+1', rendered)
        self.assertNotIn('${prefix}', rendered); self.assertIn('fixture.boolean.txt >> fixture.saf', rendered)
        with self.assertRaises(probe.InfrastructureError): probe.extract_converter(actual + actual)

    def test_missing_boundary_read_or_changed_coordinate_in_BAM_rejects(self):
        records = '\n'.join(line for line in probe.SAM.splitlines() if not line.startswith('@')) + '\n'
        self.assertTrue(probe.exact_bam_records(records)['passed'])
        self.assertFalse(probe.exact_bam_records('\n'.join(records.splitlines()[:-1]))['passed'])
        self.assertFalse(probe.exact_bam_records(records.replace('\t100\t60\t', '\t101\t60\t'))['passed'])
        self.assertFalse(probe.exact_bam_records(records + records.splitlines(True)[0])['passed'])

    def test_corrupt_successful_count_output_is_scientific_failure(self):
        valid = 'Geneid\tChr\tStart\tEnd\tStrand\tLength\ttiny.bam\nInterval_zero\tchrPoc\t1\t10\t+\t10\t1\n'
        for text in ['', 'wrong header\n', valid.replace('\t1\n', '\tmissing\n'), valid + valid.splitlines(True)[-1]]:
            with self.assertRaises(probe.ScientificFailure): probe.read_counts(text)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'counts.tsv'
            with self.assertRaises(probe.ScientificFailure): probe.produced_text(path)
            path.write_bytes(b'\xff')
            with self.assertRaises(probe.ScientificFailure): probe.produced_text(path)
        with self.assertRaises(probe.InfrastructureError): probe.extract_converter('unsupported production template')

    def test_real_count_contract_requires_boundary_exclusion_and_positive_controls(self):
        self.assertEqual(probe.count_checks({'Interval_zero': 1, 'Interval_nonzero': 1})['status'], 'PASS')
        for counts in [{'Interval_zero': 1, 'Interval_nonzero': 2}, {'Interval_zero': 0, 'Interval_nonzero': 1}, {'Interval_nonzero': 1}]:
            self.assertEqual(probe.count_checks(counts)['status'], 'FAIL')
        text = '# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ttiny.bam\nInterval_zero\tchrPoc\t1\t10\t+\t10\t1\n'
        self.assertEqual(probe.read_counts(text), {'Interval_zero': 1})
        with self.assertRaises(probe.ScientificFailure): probe.read_counts(text + text.splitlines(True)[-1])


if __name__ == '__main__': unittest.main()
