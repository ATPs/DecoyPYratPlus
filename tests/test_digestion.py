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
    def run_script_cli(self, *args):
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "DecoyPYratPlus" / "digestion.py"), *args],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout

    def run_cli(self, *args):
        result = subprocess.run(
            [sys.executable, "-m", "DecoyPYratPlus.digestion", *args],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout

    def test_help_contains_examples_and_enzyme_presets(self):
        help_output = self.run_script_cli("--help")
        self.assertIn("keep I and L distinct", help_output)
        self.assertIn("python DecoyPYratPlus/digestion.py --sequence AKRPQK -l 2", help_output)
        self.assertIn("--header-template \"{protein_id}|pep{index}\"", help_output)
        self.assertIn("Cut_everywhere", help_output)
        self.assertIn("No_cut", help_output)
        self.assertIn("Comet-style enzyme presets", help_output)
        self.assertIn(">sp|P12345\\tAK", help_output)

    def test_default_sequence_cli_output_formats(self):
        tsv_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "--enzyme",
            "Trypsin",
            "-l",
            "2",
            "--output-format",
            "tsv",
        )
        self.assertEqual(tsv_output, "sequence\tAK\nsequence\tRPQK\n")

        peptide_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "--enzyme",
            "Trypsin",
            "-l",
            "2",
            "--output-format",
            "peptide",
        )
        self.assertEqual(peptide_output, "AK\nRPQK\n")

        fasta_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "--enzyme",
            "Trypsin",
            "-l",
            "2",
            "--output-format",
            "fasta",
        )
        self.assertEqual(fasta_output, ">sequence_1\nAK\n>sequence_2\nRPQK\n")

    def test_isobaric_default_and_override(self):
        peptide_output = self.run_cli(
            "--sequence",
            "AILK",
            "--enzyme",
            "No_cut",
            "-l",
            "1",
            "--output-format",
            "peptide",
        )
        self.assertEqual(peptide_output, "AILK\n")

        peptide_output_isobaric = self.run_cli(
            "--sequence",
            "AILK",
            "--enzyme",
            "No_cut",
            "--isobaric",
            "-l",
            "1",
            "--output-format",
            "peptide",
        )
        self.assertEqual(peptide_output_isobaric, "ALLK\n")

    def test_enzyme_presets_and_manual_override(self):
        no_cut_output = self.run_cli(
            "--sequence",
            "AILK",
            "--enzyme",
            "No_cut",
            "-l",
            "1",
            "--output-format",
            "peptide",
        )
        self.assertEqual(no_cut_output, "AILK\n")

        cut_everywhere_output = self.run_cli(
            "--sequence",
            "AKRPQ",
            "--enzyme",
            "Cut_everywhere",
            "-l",
            "2",
            "-M",
            "3",
            "--output-format",
            "peptide",
        )
        self.assertEqual(cut_everywhere_output, "AK\nKR\nRP\nPQ\nAKR\nKRP\nRPQ\n")

        trypsin_output = self.run_cli(
            "--sequence",
            "AKRPQK",
            "--enzyme",
            "Trypsin",
            "-l",
            "2",
            "--output-format",
            "peptide",
        )
        self.assertEqual(trypsin_output, "AK\nRPQK\n")

        lys_n_output = self.run_cli(
            "--sequence",
            "KAAKBB",
            "--enzyme",
            "Lys_N",
            "-l",
            "1",
            "--output-format",
            "peptide",
        )
        self.assertEqual(lys_n_output, "KAA\nKBB\n")

        override_output = self.run_cli(
            "--sequence",
            "KAAKBB",
            "--enzyme",
            "Cut_everywhere",
            "-c",
            "K",
            "-a",
            "",
            "-p",
            "c",
            "-l",
            "1",
            "--output-format",
            "peptide",
        )
        self.assertEqual(override_output, "K\nAAK\nBB\n")

    def test_fasta_cli_output_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "proteins.fa"
            fasta_path.write_text(">prot1\nAKBKCK\n", encoding="ascii")

            tsv_output = self.run_cli(
                str(fasta_path),
                "--method",
                "trypsin",
                "--enzyme",
                "Trypsin/P",
                "-c",
                "K",
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
                "--enzyme",
                "Trypsin/P",
                "-c",
                "K",
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
                "--enzyme",
                "Trypsin/P",
                "-c",
                "K",
                "-L",
                "1",
                "-l",
                "2",
                "-M",
                "4",
                "--output-format",
                "fasta",
                "--header-template",
                "{protein_id}|pep{index}",
            )
            self.assertEqual(
                fasta_output,
                ">prot1|pep1\nAK\n>prot1|pep2\nBK\n>prot1|pep3\nCK\n>prot1|pep4\nAKBK\n>prot1|pep5\nBKCK\n",
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

    def test_main_decoypyratplus_default_no_isobaric_behavior_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_fasta = Path(tmpdir) / "input.fa"
            default_output = Path(tmpdir) / "default.fa"
            no_isobaric_output = Path(tmpdir) / "no_isobaric.fa"
            input_fasta.write_text(">p1\nAIIK\n", encoding="ascii")

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "DecoyPYratPlus.decoyPYratPlus",
                    str(input_fasta),
                    "-o",
                    str(default_output),
                    "--keep_names",
                    "--do_not_shuffle",
                    "--do_not_switch",
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("KLLA", default_output.read_text(encoding="ascii"))

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "DecoyPYratPlus.decoyPYratPlus",
                    str(input_fasta),
                    "-o",
                    str(no_isobaric_output),
                    "--keep_names",
                    "--do_not_shuffle",
                    "--do_not_switch",
                    "--no_isobaric",
                ],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("KIIA", no_isobaric_output.read_text(encoding="ascii"))


if __name__ == "__main__":
    unittest.main()
