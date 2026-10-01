#!/usr/bin/env bash
# Before/after table for one tier, against the snapshot from 00_snapshot_baseline.sh.
#
# Usage: ./99_compare_to_baseline.sh tier3
set -euo pipefail

TIER="${1:-}"
if [ -z "$TIER" ]; then
    echo "Usage: $0 <tier2|tier3|tier4|tier5>" >&2
    exit 1
fi

RUN=NOVA_260909_LA
OUT="/scratch/pawsey1348/$USER/${RUN}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
BASELINE="$OUT/rerun_baseline"
SAMPLESHEET="$HERE/samplesheets/${TIER}_samplesheet.csv"

for p in "$BASELINE" "$SAMPLESHEET"; do
    [ -e "$p" ] || { echo "ERROR: $p not found" >&2; exit 1; }
done

python3 - "$OUT" "$BASELINE" "$SAMPLESHEET" <<'PY'
import csv, glob, json, os, re, sys

out, baseline, ss = sys.argv[1], sys.argv[2], sys.argv[3]
samples = [r[0] for r in list(csv.reader(open(ss)))[1:] if r and r[0]]

def busco_c(d):
    for f in glob.glob(f"{d}/*short_summary.txt") + glob.glob(f"{d}/assemblies/genome/busco/*short_summary.txt"):
        m = re.search(r"C:([0-9.]+)%", open(f).read())
        if m:
            return float(m.group(1))
    return None

def seqkit(d):
    for f in glob.glob(f"{d}/*seqkit_stats.tsv") + glob.glob(f"{d}/assemblies/genome/seqkit/*seqkit_stats.tsv"):
        rows = [l.split('\t') for l in open(f).read().splitlines() if l]
        if len(rows) > 1:
            r = rows[1]
            return int(r[4].replace(',', '')), int(r[3].replace(',', '')), int(r[12].replace(',', ''))
    return None, None, None

def _summary(d):
    for f in glob.glob(f"{d}/*_summary.txt") + glob.glob(f"{d}/kmers/*genomescope/*_summary.txt"):
        if "genomescope" in os.path.basename(f) or "Model Fit" in open(f).read():
            return f
    return None

def _bound(token):
    """A GenomeScope bound. "Inf" comes back as None -- never as the other bound. The
    earlier version of this matched only digits, so a diverged fit (OG3043 reported
    "2,499,788 bp - Inf bp") silently displayed its LOWER bound as the estimate."""
    if re.fullmatch(r"-?(inf|na|nan)", token.strip(), flags=re.I):
        return None
    try:
        return float(token.replace(",", ""))
    except ValueError:
        return None

def _prop(d, label, which):
    f = _summary(d)
    if not f:
        return None
    m = re.search(rf"^{label}\s+(\S+)\s*(?:bp)?\s+(\S+)\s*(?:bp)?\s*$",
                  open(f).read(), flags=re.M)
    if not m:
        return None
    return _bound(m.group(1 if which == "min" else 2).rstrip("%"))

# The UPPER model fit is the fit to the modelled region and is what the gate thresholds.
# The lower one is a residual over the -m window: it falls whenever that ceiling rises,
# which is exactly the artifact that made a parameter change look like a model regression.
def gs_fit(d):
    return _prop(d, "Model Fit", "max")

def gs_size(d):
    return _prop(d, "Genome Haploid Length", "max")

# Unique length is the component the model actually constrains. When a size estimate moves
# but this does not, the extra basepairs came out of the histogram tail, not the genome.
def gs_unique(d):
    return _prop(d, "Genome Unique Length", "max")

def cov_status(d):
    for f in glob.glob(f"{d}/*_coverage_summary.json") + glob.glob(f"{d}/coverage/*_coverage_summary.json"):
        j = json.load(open(f))
        return j.get("coverage_status"), j.get("genome_size_flags", "")
    return None, ""

hdr = (f"{'SAMPLE':<8} {'BUSCO_C':>15} {'ASM_Mb':>17} {'N50':>15} {'CTGS':>17} "
       f"{'GS_fitmax':>15} {'GS_uniq_Mb':>17} {'GS_Mb':>17}")
print(hdr)
print("-" * len(hdr))

def cell(before, after, fmt="{:.1f}", width=15):
    b = fmt.format(before) if before is not None else "NA"
    a = fmt.format(after) if after is not None else "NA"
    return f"{b}->{a}".rjust(width)

def mb(value):
    return value / 1e6 if value else None

for s in samples:
    bd = f"{baseline}/{s}"
    ad = f"{out}/draftgenomes/{s}"
    if not os.path.isdir(bd):
        print(f"{s:<8} (no baseline)")
        continue
    b_asm, b_ctg, b_n50 = seqkit(bd)
    a_asm, a_ctg, a_n50 = seqkit(ad)
    print(f"{s:<8}"
          f"{cell(busco_c(bd), busco_c(ad))}"
          f"{cell(b_asm/1e6 if b_asm else None, a_asm/1e6 if a_asm else None, width=17)}"
          f"{cell(b_n50, a_n50, '{:.0f}')}"
          f"{cell(b_ctg, a_ctg, '{:.0f}', 17)}"
          f"{cell(gs_fit(bd), gs_fit(ad))}"
          f"{cell(mb(gs_unique(bd)), mb(gs_unique(ad)), '{:.0f}', 17)}"
          f"{cell(mb(gs_size(bd)), mb(gs_size(ad)), '{:.0f}', 17)}")

print()
print("Coverage verdicts after the re-run:")
for s in samples:
    st, flags = cov_status(f"{out}/draftgenomes/{s}")
    if st:
        print(f"  {s:<8} {st}{('  [' + flags + ']') if flags else ''}")
PY
