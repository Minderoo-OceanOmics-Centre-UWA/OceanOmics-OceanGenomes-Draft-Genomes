#!/usr/bin/env python3
"""Plan and adjudicate a GenomeScope reseed.

Two subcommands, matching the two points the module needs an answer:

  plan    should a reseed be attempted, and with what integer seed?
  choose  given both fits, which one does the pipeline keep?

The judgement itself lives in genomescope_reliability.choose_reseeded_fit so it can be
tested without running GenomeScope.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genomescope_reliability as gate


def load(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {}


def plan(args):
    """Print the integer seed, or "skip"."""
    if str(args.enabled).lower() not in ("true", "1", "yes"):
        print("skip")
        return 0

    coverage = load(args.coverage_json)
    depth = load(args.kmer_depth_json)

    # Rescue a fit that was flagged, or one whose fitted coverage already disagrees with
    # the measurement. The verdict alone is not enough: it is read from the PROVISIONAL
    # summary, computed before the assembly exists, so it cannot see either the assembly
    # cross-check or the lambda ratio. OG2617 was reported "no reseed attempted" on a
    # provisional pass while its fitted/measured ratio was 1.95.
    ratio = depth.get("kmercov_over_lambda_depth")
    ratio = max(ratio, 1 / ratio) if ratio else None
    disagrees = ratio is not None and ratio > args.max_lambda_ratio

    if coverage.get("genome_size_reliable", False) and not disagrees:
        print("skip")
        return 0

    lam = depth.get("lambda_depth")
    if not lam or lam < args.min_seedable_lambda:
        # Below this there is no host peak to find. Forcing the optimiser there produces
        # fits of a few percent with non-finite genome lengths, which is worse than the
        # flagged fit it would replace.
        print("skip")
        return 0

    # -l takes an integer only.
    print(max(1, int(round(lam))))
    return 0


def choose(args):
    original = {"summary": gate.parse_summary(args.original_summary),
                "model": gate.parse_model(args.original_model)}

    attempted = args.seed != "skip"
    reseeded = {"summary": {}, "model": {}}
    if attempted and os.path.exists(args.reseeded_summary):
        reseeded = {"summary": gate.parse_summary(args.reseeded_summary),
                    "model": gate.parse_model(args.reseeded_model)}

    lam = load(args.kmer_depth_json).get("lambda_depth")

    if not attempted:
        winner, reason = "original", "no reseed attempted"
    elif not reseeded["summary"]:
        winner, reason = "original", "reseeded run produced no summary"
    else:
        winner, reason = gate.choose_reseeded_fit(
            original, reseeded, lam, {"min_reseed_fit_gain": args.min_reseed_fit_gain})

    decision = {
        "seed": None if not attempted else int(args.seed),
        "lambda_depth": lam,
        "winner": winner,
        "reason": reason,
        "original_model_fit_full": original["summary"].get("model_fit_full"),
        "original_kmercov": original["model"].get("kmercov"),
        "original_genome_size": original["summary"].get("haploid_max"),
        "reseeded_model_fit_full": reseeded["summary"].get("model_fit_full"),
        "reseeded_kmercov": reseeded["model"].get("kmercov"),
        "reseeded_genome_size": reseeded["summary"].get("haploid_max"),
    }
    with open(args.out, "w") as handle:
        json.dump(decision, handle, indent=2)
        handle.write("\n")
    print(f"{winner}: {reason}", file=sys.stderr)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("plan")
    p.add_argument("--coverage-json", required=True)
    p.add_argument("--kmer-depth-json", required=True)
    p.add_argument("--min-seedable-lambda", type=float, default=5.0)
    p.add_argument("--max-lambda-ratio", type=float, default=1.5,
                   help="reseed when fitted/measured coverage exceeds this, even if the "
                        "provisional verdict passed")
    p.add_argument("--enabled", default="true")
    p.set_defaults(func=plan)

    c = sub.add_parser("choose")
    c.add_argument("--coverage-json", required=True)
    c.add_argument("--kmer-depth-json", required=True)
    c.add_argument("--original-summary", required=True)
    c.add_argument("--original-model", required=True)
    c.add_argument("--reseeded-summary", required=True)
    c.add_argument("--reseeded-model", required=True)
    c.add_argument("--seed", required=True)
    c.add_argument("--min-reseed-fit-gain", type=float, default=2.0)
    c.add_argument("--out", required=True)
    c.set_defaults(func=choose)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
