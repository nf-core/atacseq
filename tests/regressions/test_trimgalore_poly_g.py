import gzip
import tempfile
import unittest
from pathlib import Path

from trimgalore_poly_g import assess, expected_records, output_identity, write_reads


class TrimGaloreControls(unittest.TestCase):
    def test_original_records_pass_and_order_does_not_select_a_winner(self):
        with tempfile.TemporaryDirectory() as folder:
            for mate in (1, 2):
                path = Path(folder) / f'{mate}.fq.gz'
                write_reads(path, mate)
                self.assertTrue(assess(path, mate)['passed'])
                rows = list(reversed(expected_records(mate)))
                text = ''.join(f'@{name}\n{sequence}\n+\n{quality}\n' for name, sequence, quality in rows)
                path.write_bytes(gzip.compress(text.encode(), mtime=0))
                self.assertTrue(assess(path, mate)['passed'])

    def test_clipping_and_record_corruption_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'reads.fq.gz'
            for mate in (1, 2):
                expected = expected_records(mate)
                name, sequence, quality = expected[0]
                clipped = sequence[:-24] if mate == 1 else sequence[25:]
                mutations = [expected[:-1], expected + expected[:1],
                             [(name, clipped, 'I' * len(clipped))] + expected[1:],
                             [('wrong-id', sequence, quality)] + expected[1:],
                             [(name, sequence, 'J' * len(quality))] + expected[1:]]
                for rows in mutations:
                    with self.subTest(mate=mate, records=len(rows)):
                        text = ''.join(f'@{key}\n{seq}\n+\n{qual}\n' for key, seq, qual in rows)
                        path.write_bytes(gzip.compress(text.encode(), mtime=0))
                        self.assertFalse(assess(path, mate)['passed'])

    def test_successful_empty_malformed_and_truncated_outputs_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'reads.fq.gz'
            write_reads(path, 1)
            valid = path.read_bytes()
            for content in (gzip.compress(b'', mtime=0), gzip.compress(b'@read\nACTG\n+\nIII\n', mtime=0), valid[:-5]):
                with self.subTest(bytes=len(content)):
                    path.write_bytes(content)
                    self.assertFalse(assess(path, 1)['passed'])

    def test_missing_completed_output_stays_a_scientific_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'missing.fq.gz'
            self.assertFalse(assess(path, 1)['passed'])
            self.assertIn('read_error', output_identity(path))


if __name__ == '__main__':
    unittest.main()
