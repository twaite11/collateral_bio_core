"""
Tests for SRA Magic-BLAST miner: command build, magicblast -help options, tabular parsing.

Run from repo root:
  python -m unittest tests.test_sra_blast_miner -v

To see exactly what command would be run (and what magicblast supports):
  python -m unittest tests.test_sra_blast_miner.TestMagicBlastCommand -v
  python -m unittest tests.test_sra_blast_miner.TestMagicBlastHelp -v

On the server (with magicblast installed), run the help test to see supported options:
  RUN_MAGICBLAST_HELP=1 python -m unittest tests.test_sra_blast_miner.TestMagicBlastHelp -v
If that test fails, the message shows which of -db, -sra, -sra_batch, -outfmt are missing.
"""
import os
import subprocess
import tempfile
import unittest
import unittest.mock
from pathlib import Path

# Allow running from repo root or from tests/
_root = Path(__file__).resolve().parents[1]
if str(_root) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_root))

from modules.mining.sra_blast_miner import (
    build_magicblast_command,
    load_ref_seqs,
    parse_magicblast_tabular,
    run_magic_blast,
)


def _magicblast_available():
    try:
        subprocess.run(["magicblast", "-help"], capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


class TestMagicBlastCommand(unittest.TestCase):
    """Verify the exact command and cwd built for Magic-BLAST."""

    def test_single_accession_command(self):
        cmd, cwd, batch_file = build_magicblast_command(
            ["SRR123"],
            "/data/raw/_cas13_ref_db",
            "/data/raw/out.tsv",
            magicblast_cmd="magicblast",
            num_threads=2,
        )
        self.assertIsNone(batch_file)
        self.assertIn("magicblast", cmd[0])
        self.assertEqual(cmd[1], "-sra")
        self.assertEqual(cmd[2], "SRR123")
        idx_db = cmd.index("-db")
        self.assertEqual(cmd[idx_db + 1], "_cas13_ref_db", " -db should be basename only")
        idx_out = cmd.index("-out")
        self.assertEqual(cmd[idx_out + 1], "out.tsv", " -out should be filename only")
        idx_fmt = cmd.index("-outfmt")
        self.assertEqual(cmd[idx_fmt + 1], "tabular")
        self.assertEqual(cwd, str(Path("/data/raw/_cas13_ref_db").parent.resolve()))

    def test_batch_accession_command(self):
        cmd, cwd, batch_file = build_magicblast_command(
            ["SRR1", "SRR2"],
            "/workspace/senary_bio_core/data/raw_sequences/_cas13_ref_db",
            "/workspace/senary_bio_core/data/raw_sequences/_magicblast_batch_0.tsv",
            magicblast_cmd="/opt/magicblast/bin/magicblast",
            num_threads=4,
        )
        self.assertIsNotNone(batch_file)
        self.assertTrue(batch_file.endswith(".txt"))
        self.assertIn("_sra_batch_", batch_file)
        self.assertEqual(cmd[0], "/opt/magicblast/bin/magicblast")
        self.assertEqual(cmd[1], "-sra_batch")
        # -sra_batch value must be filename only (so it's found in cwd)
        self.assertEqual(cmd[2], "_sra_batch__magicblast_batch_0.txt")
        idx_db = cmd.index("-db")
        self.assertEqual(cmd[idx_db + 1], "_cas13_ref_db")
        self.assertNotIn("/", cmd[idx_db + 1], " -db must be basename, no path")
        self.assertEqual(Path(cwd).name, "raw_sequences", "cwd should be the DB parent directory")

    def test_run_magic_blast_writes_batch_and_uses_cwd(self):
        """With mocked subprocess, verify run_magic_blast builds correct cmd and cwd."""
        tmp = _root / "_test_magicblast_tmp"
        tmp.mkdir(parents=True, exist_ok=True)
        try:
            db_path = str(tmp / "_cas13_ref_db")
            out_tsv = str(tmp / "batch_0.tsv")

            with unittest.mock.patch("modules.mining.sra_blast_miner.subprocess.run") as m_run:
                m_run.return_value = unittest.mock.Mock(returncode=0, stderr="", stdout="")
                ok = run_magic_blast(
                    ["SRR_A", "SRR_B"],
                    db_path,
                    out_tsv,
                    magicblast_cmd="magicblast",
                    num_threads=1,
                )
            self.assertTrue(ok)
            m_run.assert_called_once()
            call_args = m_run.call_args[0][0]
            call_kw = m_run.call_args[1]
            self.assertEqual(call_args[0], "magicblast")
            self.assertEqual(call_args[1], "-sra_batch")
            batch_name = call_args[2]
            self.assertTrue(batch_name.endswith(".txt"))
            self.assertEqual(call_args[call_args.index("-db") + 1], "_cas13_ref_db")
            self.assertEqual(call_args[call_args.index("-out") + 1], "batch_0.tsv")
            self.assertEqual(Path(call_kw.get("cwd")).resolve(), tmp.resolve())
            # Batch filename in cmd must be the one we write (miner deletes it in finally)
            self.assertEqual(batch_name, "_sra_batch_batch_0.txt")
        finally:
            for f in tmp.glob("*"):
                try:
                    f.unlink()
                except OSError:
                    pass
            try:
                tmp.rmdir()
            except OSError:
                pass


class TestMagicBlastHelp(unittest.TestCase):
    """Run magicblast -help (or -h) if available and check which options are supported."""

    @unittest.skipUnless(
        os.environ.get("RUN_MAGICBLAST_HELP") or _magicblast_available(),
        "Set RUN_MAGICBLAST_HELP=1 or have magicblast on PATH to run",
    )
    def test_magicblast_help_lists_required_options(self):
        """Fail with a clear message if -db, -sra, -outfmt are missing from help."""
        result = subprocess.run(
            ["magicblast", "-help"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        help_text = (result.stdout or "") + (result.stderr or "")
        required = ["-db", "-out"]
        sra_options = ["-sra", "-sra_batch"]
        format_related = ["-outfmt", "tabular", "sam"]
        missing = [r for r in required if r not in help_text]
        missing_sra = [s for s in sra_options if s not in help_text]
        any_fmt = any(f in help_text for f in format_related)
        if missing:
            self.fail(f"magicblast -help does not show required options: {missing}\nHelp snippet:\n{help_text[:1500]}")
        if missing_sra:
            self.fail(
                f"magicblast -help does not show SRA options {missing_sra}. "
                f"Your build may not support SRA mode.\nHelp snippet:\n{help_text[:1500]}"
            )
        if not any_fmt:
            self.fail(
                "magicblast -help does not show -outfmt / tabular / sam. "
                f"Output format may differ.\nHelp snippet:\n{help_text[:1500]}"
            )
        # Optional: print what we found so user can compare
        print(" [magicblast -help] Found: -db, -out, SRA options, and format-related text.")

    def test_magicblast_command_echo(self):
        """Print the exact command that would be run (no subprocess). Use for manual debugging."""
        cmd, cwd, batch_file = build_magicblast_command(
            ["SRR1237994", "SRR1237993"],
            "/workspace/senary_bio_core/data/raw_sequences/_cas13_ref_db",
            "/workspace/senary_bio_core/data/raw_sequences/_magicblast_batch_0.tsv",
            magicblast_cmd="magicblast",
            num_threads=4,
        )
        # So user can copy-paste and run on their server
        cmd_str = " ".join(cmd)
        print(f"\n[Echo] cwd={cwd}")
        print(f"[Echo] cmd={cmd_str}")
        if batch_file:
            print(f"[Echo] batch_file (create before run): {batch_file}")
        self.assertIn("-db", cmd)
        self.assertIn("-outfmt", cmd)


class TestParseMagicBlastTabular(unittest.TestCase):
    """Test tabular parsing with sample data."""

    def test_parse_tabular_forward_hit(self):
        ref_seqs = {"ref1": "ATGAAACCCGGGTTTTAAACCCGGGTTTTAAATAA"}  # 33 nt, simple
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tab", delete=False) as f:
            # Magic-BLAST tabular: qseqid, sseqid, pident, 4-6 unused, 7 qstart, 8 qend, 9 sstart, 10 send, ...
            f.write("read1\tref1\t95\t0\t0\t0\t1\t30\t1\t30\t0\t0\t100\t+\t+\t30\t30\n")
            f.flush()
            rows = parse_magicblast_tabular(f.name, ref_seqs)
        os.unlink(f.name)
        self.assertEqual(len(rows), 1)
        qseqid, sseqid, protein, sseq_nt, sstart, send = rows[0]
        self.assertEqual(qseqid, "read1")
        self.assertEqual(sseqid, "ref1")
        self.assertEqual(sstart, 1)
        self.assertEqual(send, 30)
        self.assertEqual(len(sseq_nt), 30)
        # Forward: ATG (M) then 27 nt = 9 aa
        self.assertTrue(protein.startswith("M") or len(protein) >= 9)

    def test_parse_tabular_skips_empty_and_comments(self):
        ref_seqs = {"r": "ATG" * 10}
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tab", delete=False) as f:
            f.write("# comment\n")
            f.write("q\tr\t90\t0\t0\t0\t1\t9\t1\t9\t0\t0\t50\t+\t+\t9\t9\n")
            f.write("\n")
            f.flush()
            rows = parse_magicblast_tabular(f.name, ref_seqs)
        os.unlink(f.name)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], "q")
