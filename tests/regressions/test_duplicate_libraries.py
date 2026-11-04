"""Exercise semantic membership failures independently of snapshot contents."""
import importlib.util
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('duplicate_libraries', Path(__file__).with_name('duplicate_libraries.py'))
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class DuplicateMembershipControls(unittest.TestCase):
    def fixture(self):
        groups = {run: {'--RGLB': 'LIB_A_REP1' if i < 2 else 'LIB_B_REP1', '--RGSM': 'LIB_A_REP1' if i < 2 else 'LIB_B_REP1'}
                  for i, run in enumerate(guard.RUNS)}
        header = {run: {'ID': run, 'LB': group['--RGLB'], 'SM': group['--RGSM']} for run, group in groups.items()}
        names = list(guard.NAMES.values())
        records = [{'name': name, 'mate': mate, 'duplicate': name == guard.NAMES['a_duplicate'], 'rg': guard.READ_RUN[name],
                    'reference': 'chrPoc', 'cigar': '20M', 'position': guard.READ_START[name] + (50 if mate == 2 else 0),
                    'flag': (99 if mate == 1 else 147) | (1024 if name == guard.NAMES['a_duplicate'] else 0)}
                   for name in names for mate in (1, 2)]
        return header, records

    def test_exact_membership_passes(self):
        self.assertEqual(guard.scientific_checks(*self.fixture())['status'], 'PASS')

    def test_missing_technical_duplicate_is_a_failure(self):
        header, records = self.fixture()
        for record in records:
            record['duplicate'] = False
        self.assertEqual(guard.scientific_checks(header, records)['status'], 'FAIL')

    def test_independent_library_overdeduplication_is_a_failure(self):
        header, records = self.fixture()
        for record in records:
            if record['name'] == guard.NAMES['b_independent']:
                record['duplicate'] = True
        result = guard.scientific_checks(header, records)
        self.assertFalse(result['checks']['independent_biological_library_same_coordinates_retained'])
        self.assertEqual(result['status'], 'FAIL')

    def test_dropped_pair_cannot_masquerade_as_duplicate_removal(self):
        header, records = self.fixture()
        records = [r for r in records if r['name'] != guard.NAMES['a_duplicate']]
        self.assertEqual(guard.scientific_checks(header, records)['status'], 'FAIL')

    def test_single_mate_flag_is_not_a_duplicate_pair(self):
        header, records = self.fixture()
        next(r for r in records if r['name'] == guard.NAMES['a_duplicate'])['duplicate'] = False
        self.assertEqual(guard.scientific_checks(header, records)['status'], 'FAIL')

    def test_parser_preserves_duplicate_bit_and_mates(self):
        sam = '@RG\tID:r1\tLB:l1\tSM:s1\nq\t1123\tchrPoc\t100\t60\t20M\t=\t150\t70\tAAAA\tIIII\tRG:Z:r1\n'
        header, rows = guard.parse_sam(sam)
        self.assertEqual(header['r1']['LB'], 'l1')
        self.assertTrue(rows[0]['duplicate'])
        self.assertEqual(rows[0]['mate'], 1)

    def test_malformed_produced_alignment_is_scientific_failure(self):
        with self.assertRaises(guard.ScientificFailure):
            guard.parse_sam('read\tbad-flag\tchrPoc\t100\t60\t20M\t=\t150\t70\tAAAA\tIIII\n')

    def test_duplicate_produced_read_group_is_scientific_failure(self):
        with self.assertRaises(guard.ScientificFailure):
            guard.parse_sam('@RG\tID:r1\tLB:l1\n@RG\tID:r1\tLB:l2\n')


if __name__ == '__main__':
    unittest.main()
