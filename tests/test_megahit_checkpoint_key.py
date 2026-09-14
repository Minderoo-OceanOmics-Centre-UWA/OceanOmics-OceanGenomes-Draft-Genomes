"""Tests for MEGAHIT checkpoint fingerprinting.

Two layers:

* the fingerprint itself (bin/megahit_checkpoint_key.sh), which decides whether an
  existing checkpoint belongs to the current inputs;
* the decision logic in modules/local/megahit/main.nf, exercised by rendering its
  script block into plain bash and running it against a stand-in megahit.

The second layer matters more than it looks. The checkpoint directory lives outside the
work directory so a killed assembly can resume, which puts it beyond Nextflow's caching:
nothing else in the pipeline will notice if a stale assembly is republished.
"""

import os
import re
import shutil
import subprocess
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KEY_SCRIPT = ROOT / "bin" / "megahit_checkpoint_key.sh"
MODULE = ROOT / "modules" / "local" / "megahit" / "main.nf"

FAKE_MEGAHIT = """#!/bin/bash
# Stand-in for megahit: enough of its behaviour to exercise the checkpoint logic.
if [ "$1" = "-v" ]; then echo "MEGAHIT v1.2.9"; exit 0; fi
out=""; prefix="out"; cont=false
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out="$2"; shift 2;;
    --out-prefix) prefix="$2"; shift 2;;
    --continue) cont=true; shift;;
    *) shift;;
  esac
done
if $cont; then
  [ -f "$out/opts.txt" ] && prefix=$(cat "$out/opts.txt")
  [ -n "${FAKE_CONTINUE_FAIL:-}" ] && { echo "continue failed" >&2; exit 1; }
else
  # real megahit refuses to write into an -o that already exists
  if [ -d "$out" ]; then echo "ERROR: $out already exists" >&2; exit 1; fi
  mkdir -p "$out"
  echo checkpoint > "$out/checkpoints.txt"
  echo "$prefix" > "$out/opts.txt"
  [ -n "${FAKE_FRESH_FAIL:-}" ] && { echo "crashed" >&2; exit 1; }
fi
mkdir -p "$out/intermediate_contigs" "$out/tmp"
echo graph > "$out/intermediate_contigs/k21.contigs.fa"
echo scratch > "$out/tmp/scratch"
printf '>contig_%s\\nACGT\\n' "${FAKE_TAG:-A}" > "$out/$prefix.contigs.fa"
touch "$out/done"
"""


def key(*reads, version="1.2.9", args="", prefix="OG1", cwd=None, manifest=False):
    cmd = [str(KEY_SCRIPT), "--version", version, "--args", args, "--prefix", prefix]
    if manifest:
        cmd.append("--manifest")
    cmd += list(reads)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)


class TestCheckpointKey(unittest.TestCase):
    """The fingerprint keys what changes the assembly, and nothing else."""

    def setUp(self):
        self.dir = Path(os.environ.get("TMPDIR", "/tmp")) / f"mhkey_{os.getpid()}_{id(self)}"
        self.dir.mkdir(parents=True, exist_ok=True)
        (self.dir / "OG1.fastp.R1.fastq.gz").write_text("aaaaaaaa")
        (self.dir / "OG1.fastp.R2.fastq.gz").write_text("bbbbbbbb")
        (self.dir / "OG1.kraken2filt.R1.fastq.gz").write_text("aaaa")
        (self.dir / "OG1.kraken2filt.R2.fastq.gz").write_text("bbbb")
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)

    def fastp(self):
        return ["OG1.fastp.R1.fastq.gz", "OG1.fastp.R2.fastq.gz"]

    def k(self, *reads, **kw):
        kw.setdefault("cwd", self.dir)
        res = key(*reads, **kw)
        self.assertEqual(res.returncode, 0, res.stderr)
        return res.stdout.strip()

    def test_stable_across_channel_order(self):
        a, b = self.fastp()
        self.assertEqual(self.k(a, b), self.k(b, a))

    def test_resource_flags_are_excluded(self):
        # -m and -t change on every retry via conf/base.config and do not change the
        # contigs. Keying on them would wipe the checkpoint on the retry it exists for.
        base = self.k(*self.fastp())
        self.assertEqual(base, self.k(*self.fastp(), args="-m 100000 -t 16"))
        self.assertEqual(base, self.k(*self.fastp(), args="-m 230000 -t 32"))
        self.assertEqual(base, self.k(*self.fastp(), args="--memory=0.9 --num-cpu-threads=8"))

    def test_assembly_args_change_the_key(self):
        self.assertNotEqual(self.k(*self.fastp()), self.k(*self.fastp(), args="--k-min 27"))

    def test_different_reads_change_the_key(self):
        # The case this exists for: kraken2 filtering inserted upstream of assembly.
        self.assertNotEqual(
            self.k(*self.fastp()),
            self.k("OG1.kraken2filt.R1.fastq.gz", "OG1.kraken2filt.R2.fastq.gz"),
        )

    def test_read_content_change_changes_the_key(self):
        before = self.k(*self.fastp())
        (self.dir / "OG1.fastp.R1.fastq.gz").write_text("aaaaaaaaa")
        self.assertNotEqual(before, self.k(*self.fastp()))

    def test_version_and_prefix_change_the_key(self):
        base = self.k(*self.fastp())
        self.assertNotEqual(base, self.k(*self.fastp(), version="1.3.0"))
        self.assertNotEqual(base, self.k(*self.fastp(), prefix="OG2"))

    def test_symlinked_reads_are_followed(self):
        # Reads reach a task as symlinks into the work directory.
        (self.dir / "link.R1.gz").symlink_to(self.dir / "OG1.fastp.R1.fastq.gz")
        out = self.k("link.R1.gz", manifest=True)
        self.assertIn("read=link.R1.gz\t8", out)

    def test_missing_read_fails_without_emitting_a_key(self):
        res = key("nope.fastq.gz", cwd=self.dir)
        self.assertEqual(res.returncode, 1)
        # A partial manifest would hash to a perfectly plausible key.
        self.assertEqual(res.stdout.strip(), "")
        self.assertIn("does not exist", res.stderr)


def render_script_block(**overrides):
    """Turn the module's script block into runnable bash, as Nextflow would."""
    src = MODULE.read_text()
    start = src.index('"""', src.index("    script:")) + 3
    body = textwrap.dedent(src[start:src.index('"""', start)])

    subs = {
        "${memory}": "100000000000",
        "${args}": "",
        "${prefix}": "OG1",
        "${reads_command}": "-1 OG1.R1.fq.gz -2 OG1.R2.fq.gz",
        "${read_files}": "OG1.R1.fq.gz OG1.R2.fq.gz",
        "${checkpoint_base}": "CKBASE",
        "${output_dir}": "CKBASE/OG1_megahit_out",
        "${unkeyed_policy}": "invalidate",
        "${stale_policy}": "rerun",
        "${cleanup}": "true",
        "${effective_args}": "-m 100000000000 -t 16",
        "${task.cpus}": "16",
        "${task.process}": "TEST:MEGAHIT",
    }
    subs.update(overrides)
    for token, value in subs.items():
        body = body.replace(token, value)

    # Groovy unescapes \$ and \\ in a triple-quoted string; \t and \n are left for
    # the shell.
    body = re.sub(r"\\(.)", lambda m: m.group(1) if m.group(1) in "$\\" else "\\" + m.group(1),
                  body, flags=re.S)
    unresolved = re.findall(r"\$\{(?:task|params)\.[A-Za-z_.]+\}", body)
    assert not unresolved, f"render_script_block needs substitutions for {set(unresolved)}"
    return body


class TestCheckpointDecisions(unittest.TestCase):
    """Which path the module takes for each state a checkpoint can be in."""

    def setUp(self):
        self.dir = Path(os.environ.get("TMPDIR", "/tmp")) / f"mhrun_{os.getpid()}_{id(self)}"
        (self.dir / "bin").mkdir(parents=True, exist_ok=True)
        fake = self.dir / "bin" / "megahit"
        fake.write_text(FAKE_MEGAHIT)
        fake.chmod(0o755)
        shutil.copy(KEY_SCRIPT, self.dir / "bin" / KEY_SCRIPT.name)
        (self.dir / "OG1.R1.fq.gz").write_text("aaaaaaaa")
        (self.dir / "OG1.R2.fq.gz").write_text("bbbbbbbb")
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)

    def run_task(self, expect_success=True, **kwargs):
        script = self.dir / "task.sh"
        script.write_text(render_script_block(**kwargs.pop("overrides", {})))
        env = dict(os.environ)
        env["PATH"] = f"{self.dir / 'bin'}:{env['PATH']}"
        env.update(kwargs.pop("env", {}))
        res = subprocess.run(["bash", "-ue", str(script)], cwd=self.dir, env=env,
                             capture_output=True, text=True)
        if expect_success:
            self.assertEqual(res.returncode, 0, res.stdout + res.stderr)
        return res

    @property
    def status(self):
        text = (self.dir / "OG1.megahit_checkpoint.txt").read_text()
        return dict(line.split("\t", 1) for line in text.strip().split("\n"))

    def contig(self):
        return (self.dir / "OG1.v129mh.fasta").read_text().splitlines()[0]

    def legacy_checkpoint(self, tag="OLD"):
        """A completed checkpoint with no fingerprint, as every pre-existing one is."""
        out = self.dir / "CKBASE" / "OG1_megahit_out"
        (out / "intermediate_contigs").mkdir(parents=True)
        (out / "OG1.contigs.fa").write_text(f">contig_{tag}\nACGT\n")
        (out / "checkpoints.txt").write_text("checkpoint\n")
        (out / "opts.txt").write_text("OG1\n")
        (out / "done").touch()
        return out

    def test_fresh_when_nothing_on_disk(self):
        self.run_task()
        self.assertEqual(self.status["checkpoint_decision"], "fresh")
        self.assertTrue((self.dir / "CKBASE" / "OG1_megahit_out.key").is_file())

    def test_skip_when_inputs_unchanged(self):
        self.run_task()
        first = self.status["checkpoint_key"]
        self.run_task()
        self.assertEqual(self.status["checkpoint_decision"], "skip")
        self.assertEqual(self.status["checkpoint_key"], first)

    def test_changed_reads_force_a_rebuild(self):
        self.run_task()
        (self.dir / "OG1.R1.fq.gz").write_text("aaaa")
        self.run_task(env={"FAKE_TAG": "NEW"})
        self.assertEqual(self.status["checkpoint_decision"], "invalidated-stale-fresh")
        self.assertEqual(self.contig(), ">contig_NEW")

    def test_unkeyed_checkpoint_is_discarded_by_default(self):
        self.legacy_checkpoint()
        self.run_task(env={"FAKE_TAG": "REBUILT"})
        self.assertEqual(self.status["checkpoint_decision"], "invalidated-unkeyed-fresh")
        self.assertEqual(self.contig(), ">contig_REBUILT")

    def test_unkeyed_checkpoint_can_be_adopted(self):
        self.legacy_checkpoint()
        self.run_task(overrides={"${unkeyed_policy}": "adopt"}, env={"FAKE_TAG": "REBUILT"})
        self.assertEqual(self.status["checkpoint_decision"], "adopted-unkeyed-skip")
        self.assertEqual(self.contig(), ">contig_OLD")
        self.assertTrue((self.dir / "CKBASE" / "OG1_megahit_out.key").is_file())

    def test_stale_checkpoint_can_be_archived(self):
        self.run_task()
        (self.dir / "OG1.R1.fq.gz").write_text("aaaa")
        self.run_task(overrides={"${stale_policy}": "archive"})
        archived = list((self.dir / "CKBASE").glob("OG1_megahit_out.stale.*"))
        self.assertEqual(len([p for p in archived if p.is_dir()]), 1)
        self.assertEqual(len([p for p in archived if p.suffix == ".key"]), 1)

    def test_stale_checkpoint_can_fail_the_task(self):
        self.run_task()
        (self.dir / "OG1.R1.fq.gz").write_text("aaaa")
        res = self.run_task(overrides={"${stale_policy}": "fail"}, expect_success=False)
        self.assertEqual(res.returncode, 1)
        self.assertIn("built from different inputs", res.stderr)

    def test_crashed_assembly_resumes_on_retry(self):
        # The whole point of the checkpoint: an OOM must not cost the work done so far.
        res = self.run_task(env={"FAKE_FRESH_FAIL": "1"}, expect_success=False)
        self.assertEqual(res.returncode, 1)
        self.assertTrue((self.dir / "CKBASE" / "OG1_megahit_out.key").is_file())
        self.run_task()
        self.assertEqual(self.status["checkpoint_decision"], "continue")

    def test_retry_at_higher_resources_still_resumes(self):
        self.run_task(env={"FAKE_FRESH_FAIL": "1"}, expect_success=False)
        self.run_task(overrides={"${memory}": "230000000000", "${task.cpus}": "32"})
        self.assertEqual(self.status["checkpoint_decision"], "continue")

    def test_failed_resume_falls_back_to_a_rebuild(self):
        self.run_task(env={"FAKE_FRESH_FAIL": "1"}, expect_success=False)
        self.run_task(env={"FAKE_CONTINUE_FAIL": "1", "FAKE_TAG": "REBUILT"})
        self.assertEqual(self.status["checkpoint_decision"], "resume-failed-rebuilt")
        self.assertEqual(self.contig(), ">contig_REBUILT")

    def test_done_without_contigs_is_treated_as_partial(self):
        out = self.legacy_checkpoint()
        (out / "OG1.contigs.fa").write_text("")
        self.run_task(overrides={"${unkeyed_policy}": "adopt"}, env={"FAKE_TAG": "REBUILT"})
        self.assertEqual(self.status["checkpoint_decision"], "adopted-unkeyed-continue")
        self.assertEqual(self.contig(), ">contig_REBUILT")

    def test_completed_checkpoint_is_pruned(self):
        self.run_task()
        out = self.dir / "CKBASE" / "OG1_megahit_out"
        self.assertEqual(self.status["checkpoint_pruned"], "yes")
        self.assertFalse((out / "intermediate_contigs").exists())
        self.assertTrue((out / "OG1.contigs.fa").exists())
        self.assertTrue((out / "done").exists())

    def test_pruning_can_be_disabled(self):
        self.run_task(overrides={"${cleanup}": "false"})
        self.assertEqual(self.status["checkpoint_pruned"], "no")
        self.assertTrue((self.dir / "CKBASE" / "OG1_megahit_out" / "intermediate_contigs").exists())

    def test_unknown_policies_are_rejected_before_any_work(self):
        # Checked up front, not only when a stale checkpoint happens to turn up: the
        # schema treats a bad value as a warning, so this is the only real guard.
        for override in ({"${unkeyed_policy}": "bogus"}, {"${stale_policy}": "bogus"}):
            res = self.run_task(overrides=override, expect_success=False)
            self.assertEqual(res.returncode, 1)
            self.assertIn("ERROR: unknown", res.stderr)
            self.assertFalse((self.dir / "CKBASE").exists())

    def test_status_file_is_tab_separated(self):
        self.run_task()
        text = (self.dir / "OG1.megahit_checkpoint.txt").read_text()
        self.assertIn("checkpoint_decision\tfresh", text)
        self.assertEqual(self.status["megahit_version"], "1.2.9")
        self.assertEqual(self.status["sample"], "OG1")


if __name__ == "__main__":
    unittest.main()
