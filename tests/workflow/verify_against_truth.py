#!/usr/bin/env python3
"""Check workflow output against the fixture ground truth."""
import json, sys, os
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

print()
if failures:
    print(f"FAILED ({len(failures)}):")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print(f"ALL CHECKS PASSED across {merged['sample'].nunique()} input types "
      f"({len(merged)} rows).")
