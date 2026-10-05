#!/usr/bin/env python3
"""Check workflow output against the fixture ground truth."""
import json, sys, os, shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
truth = json.load(open(os.path.join(HERE, "fixtures", "expected.json")))
expected = pd.DataFrame(truth["expected"])
merged = pd.read_table(os.path.join(HERE, "results", "scissors_counts.tsv.gz"))

print(f"reads per distinct UMI in fixtures: {truth['reads_per_umi']}")
print("(so any pipeline counting READS instead of UMIs reports 3x these numbers)\n")

failures = []
for sample in sorted(merged["sample"].unique()):
    got = merged[merged["sample"] == sample]
    cmp = expected.merge(got, on=["CBC", "ref_name"], suffixes=("_exp", "_got"), how="outer",
                         indicator=True)
    print(f"=== {sample} ===")
    print(f"  {'cell':12s} {'template':12s} {'Pos exp/got':>14s} {'Neg exp/got':>12s} "
          f"{'Rep_Index exp/got':>22s}")
    for _, r in cmp.iterrows():
        if r["_merge"] != "both":
            failures.append(f"{sample}: {r['CBC']}/{r['ref_name']} only in {r['_merge']}")
            print(f"  MISSING: {r['CBC']}/{r['ref_name']} ({r['_merge']})")
            continue
        pos_ok = int(r["Pos_exp"]) == int(r["Pos_got"])
        neg_ok = int(r["Neg_exp"]) == int(r["Neg_got"])
        ri_exp, ri_got = r["Rep_Index_exp"], r["Rep_Index_got"]
        ri_ok = (pd.isna(ri_exp) and pd.isna(ri_got)) or abs(ri_exp - ri_got) < 1e-9
        mark = "ok " if (pos_ok and neg_ok and ri_ok) else "BAD"
        print(f"  {mark} {str(r['CBC'])[:10]:10s} {r['ref_name']:12s} "
              f"{int(r['Pos_exp']):6d}/{int(r['Pos_got']):<6d} "
              f"{int(r['Neg_exp']):5d}/{int(r['Neg_got']):<5d} "
              f"{ri_exp:10.6f}/{ri_got:<10.6f}")
        if not (pos_ok and neg_ok and ri_ok):
            failures.append(f"{sample}: {r['CBC']}/{r['ref_name']} mismatch")
    print()

# The regression the original had: two templates in one cell must differ.
print("=== REGRESSION CHECK: per-template resolution within one cell ===")
print("The original pivot_table(index='CBC') gave both templates the same")
print("Rep_Index. They must differ here.\n")
for sample in sorted(merged["sample"].unique()):
    cell = merged[(merged["sample"] == sample) & (merged["CBC"].str.startswith("AAAA"))]
    reporters = cell[cell["ref_name"].isin(["eGFP", "mRuby3"])].sort_values("ref_name")
    if len(reporters) < 2:
        continue
    vals = reporters.set_index("ref_name")["Rep_Index"]
    distinct = abs(vals["eGFP"] - vals["mRuby3"]) > 1e-9
    print(f"  {sample}: eGFP={vals['eGFP']:.6f}  mRuby3={vals['mRuby3']:.6f}  "
          f"-> {'DISTINCT (correct)' if distinct else 'IDENTICAL (bug present)'}")
    if not distinct:
        failures.append(f"{sample}: per-template Rep_Index collapsed")

# End-to-end: does a slope fitted from workflow output recover the rate the
# reads were built from?
print("=== END-TO-END: fitted slope vs the rate the reads were generated at ===")
rates = truth.get("population_rates", {})
if rates and shutil.which("Rscript") is None:
    # The counts above are the workflow's own output and are fully checked
    # without R. This last check runs the R analysis on them, so it needs
    # Rscript; skip it like the runner skips missing tools, unless the caller
    # asked for missing dependencies to be failures.
    if os.environ.get("SCISSORS_REQUIRE_DEPS") == "1":
        failures.append("Rscript not on PATH; cannot run the slope check")
    else:
        print("  SKIP: Rscript not on PATH, so analysis/scissors_replication.R")
        print("  cannot be run. Install the R side of workflow/envs/scissors.yaml.")
elif rates:
    import subprocess, tempfile, csv
    figdir = tempfile.mkdtemp(prefix="scissors_fit_")
    proc = subprocess.run(
        ["Rscript", "analysis/scissors_replication.R",
         f"--counts={os.path.join(HERE, 'results', 'scissors_counts.tsv.gz')}",
         f"--figures={figdir}", "--min-umis=50"],
        capture_output=True, text=True, cwd=os.path.join(HERE, "..", ".."))
    slopes_csv = os.path.join(figdir, "replication_slopes.csv")
    if not os.path.exists(slopes_csv):
        print("  could not run analysis/scissors_replication.R:")
        print("  " + (proc.stderr or proc.stdout).strip()[-500:])
        failures.append("analysis script did not produce replication_slopes.csv")
    else:
        with open(slopes_csv) as fh:
            rows = [r for r in csv.DictReader(fh) if r["status"] == "ok"]
        print(f"  {'sample':12s} {'template':12s} {'cells':>6s} {'true':>7s} "
              f"{'fitted':>9s} {'95% CI':>22s}")
        for r in sorted(rows, key=lambda x: (x["ref_name"], x["sample"])):
            ref = r["ref_name"]
            if ref not in rates:
                continue
            true, got = rates[ref], float(r["slope"])
            lo, hi = float(r["conf_low"]), float(r["conf_high"])
            inside = lo <= true <= hi
            print(f"  {r['sample']:12s} {ref:12s} {r['n_cells']:>6s} {true:>7.3f} "
                  f"{got:>9.5f} {'[%.5f, %.5f]' % (lo, hi):>22s} "
                  f"{'ok' if inside else 'TRUE RATE OUTSIDE CI'}")
            if not inside:
                failures.append(f"{r['sample']}/{ref}: true {true} outside CI [{lo},{hi}]")
        if not rows:
            failures.append("no groups were fitted successfully")
        # The figures the script is supposed to have drawn.
        for fig in ("neg_vs_pos.pdf", "replication_rate_slope.pdf",
                    "replication_rate_index.pdf", "donor_acceptor.pdf"):
            if not os.path.exists(os.path.join(figdir, fig)):
                failures.append(f"figure not produced: {fig}")
        produced = sorted(f for f in os.listdir(figdir) if f.endswith(".pdf"))
        print(f"\n  figures produced: {', '.join(produced)}")

print()
if failures:
    print(f"FAILED ({len(failures)}):")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print(f"ALL CHECKS PASSED across {merged['sample'].nunique()} input types "
      f"({len(merged)} rows).")
