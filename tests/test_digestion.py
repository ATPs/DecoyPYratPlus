import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from DecoyPYratPlus import decoyPYratPlus
from DecoyPYratPlus import digestion

REPO_ROOT = Path(__file__).resolve().parents[1]


class DigestionFunctionTests(unittest.TestCase):
    def test_digest_cleavage_and_anti_cleavage(self):
        peptides = digestion.digest("AKRPQK", sites="KR", pos="c", no="P", min_len=2)
        self.assertEqual(peptides, ["AK", "RPQK"])

    def test_digest_n_terminal_cleavage(self):
        peptides = digestion.digest("KAAKBB", sites="K", pos="n", no="", min_len=1)
        self.assertEqual(peptides, ["KAA", "KBB"])

    def test_trypsin_missed_cleavage_and_max_length(self):
        peptides = digestion.TRYPSIN(
            "AKBKCK",
            sites="K",
            pos="c",
            no="",
            miss_cleavage=1,
            peplen_min=2,
            peplen_max=4,
        )
        self.assertEqual(peptides, ["AK", "BK", "CK", "AKBK", "BKCK"])

        peptides_short = digestion.TRYPSIN(
            "AKBKCK",
            sites="K",
            pos="c",
            no="",
            miss_cleavage=1,
            peplen_min=2,
            peplen_max=3,
        )
        self.assertEqual(peptides_short, ["AK", "BK", "CK"])

    def test_set_fast_digest_falls_back(self):
        digestion.set_fast_digest(False)
        with redirect_stdout(io.StringIO()):
            digestion.set_fast_digest(True)
        self.assertEqual(
            digestion.digest("AKRPQK", sites="KR", pos="c", no="P", min_len=2),
            ["AK", "RPQK"],
        )

    def test_decoypyratplus_imports_migrated_symbols(self):
        self.assertIs(decoyPYratPlus.digest, digestion.digest)
        self.assertIs(decoyPYratPlus.TRYPSIN, digestion.TRYPSIN)
        self.assertIs(decoyPYratPlus.get_target_peptides, digestion.get_target_peptides)
        self.assertIs(decoyPYratPlus.get_decoy_peptides, digestion.get_decoy_peptides)


class DigestionCliTests(unittest.TestCase):
    def run_cli(self, *args):
        result = subprocess.run(
            [sys.executable, "-m", "DecoyPYratPlus.digestion", *args],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout

    def test_sequence_cli_output_formats(self):
        tsv_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "-c",
            "KR",
            "-a",
            "P",
            "-l",
            "2",
            "--output-format",
            "tsv",
        )
        self.assertEqual(tsv_output, "sequence\tAK\nsequence\tRPQK\n")

        peptide_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "-c",
            "KR",
            "-a",
            "P",
            "-l",
            "2",
            "--output-format",
            "peptide",
        )
        self.assertEqual(peptide_output, "AK\nRPQK\n")

        fasta_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "-c",
            "KR",
            "-a",
            "P",
            "-l",
            "2",
            "--output-format",
            "fasta",
        )
        self.assertEqual(fasta_output, ">sequence_1\nAK\n>sequence_2\nRPQK\n")

    def test_fasta_cli_output_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "proteins.fa"
            fasta_path.write_text(">prot1\nAKBKCK\n", encoding="ascii")

            tsv_output = self.run_cli(
                str(fasta_path),
                "--method",
                "trypsin",
                "-c",
                "K",
                "-a",
                "",
                "-L",
                "1",
                "-l",
                "2",
                "-M",
                "4",
                "--output-format",
                "tsv",
            )
            self.assertEqual(tsv_output, ">prot1\tAK\n>prot1\tBK\n>prot1\tCK\n>prot1\tAKBK\n>prot1\tBKCK\n")

            peptide_output = self.run_cli(
                str(fasta_path),
                "--method",
                "trypsin",
                "-c",
                "K",
                "-a",
                "",
                "-L",
                "1",
                "-l",
                "2",
                "-M",
                "4",
                "--output-format",
                "peptide",
            )
            self.assertEqual(peptide_output, "AK\nBK\nCK\nAKBK\nBKCK\n")

            fasta_output = self.run_cli(
                str(fasta_path),
                "--method",
                "trypsin",
                "-c",
                "K",
                "-a",
                "",
                "-L",
                "1",
                "-l",
                "2",
                "-M",
                "4",
                "--output-format",
                "fasta",
            )
            self.assertEqual(
                fasta_output,
                ">prot1_1\nAK\n>prot1_2\nBK\n>prot1_3\nCK\n>prot1_4\nAKBK\n>prot1_5\nBKCK\n",
            )

    def test_legacy_script_reuses_digestion_digest(self):
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    'import sys; '
                    'sys.argv = ["decoyPYrat.py", "dummy.fa"]; '
                    'import DecoyPYratPlus.decoyPYrat as legacy; '
                    'import DecoyPYratPlus.digestion as digestion; '
                    'print(legacy.digest is digestion.digest)'
                ),
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(result.stdout.strip(), "True")


if __name__ == "__main__":
    unittest.main()
