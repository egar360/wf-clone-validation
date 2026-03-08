#!/usr/bin/env python3
import argparse
import os
import re
import statistics
from pathlib import Path

RX = re.compile(r'_(?P<tig>[A-Z]_tig\d+?)_(?P<bp>\d+)_bp_(?P<depth>[0-9.]+)x$')

def parse_phylip(phylip: Path):
    meta = {}
    lines = phylip.read_text().splitlines()
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        nm = parts[0]
        m = RX.search(nm)
        if not m:
            continue
        tig = m.group("tig")
        bp = int(m.group("bp"))
        depth = float(m.group("depth"))
        if tig not in meta or depth > meta[tig][1]:
            meta[tig] = (bp, depth)
    return meta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--approx_size", type=int, required=True)
    ap.add_argument("--cluster_dir", default="trycycler/cluster_001")
    ap.add_argument("--phylip", default="trycycler/contigs.phylip")
    ap.add_argument("--max_keep", type=int, default=26)
    ap.add_argument("--pct_list", default="0.08,0.12,0.20")
    ap.add_argument("--min_keep", type=int, default=3)
    ap.add_argument("--reject_dir_name", default="1_contigs_rejected")
    args = ap.parse_args()

    approx_size = args.approx_size
    cluster_dir = Path(args.cluster_dir)
    contig_dir = cluster_dir / "1_contigs"
    phylip = Path(args.phylip)
    max_keep = args.max_keep
    pct_list = [float(x) for x in args.pct_list.split(",") if x.strip()]

    reject_dir = cluster_dir / args.reject_dir_name
    reject_dir.mkdir(parents=True, exist_ok=True)

    meta = parse_phylip(phylip)

    fasta_files = sorted(contig_dir.glob("*.fasta"))
    if not fasta_files:
        raise SystemExit(f"No contigs found in {contig_dir}")

    items = []
    missing = []
    for fp in fasta_files:
        tig = fp.stem
        if tig in meta:
            bp, depth = meta[tig]
            items.append((tig, fp, bp, depth))
        else:
            missing.append(tig)
            items.append((tig, fp, 0, 0.0))

    bps_nonzero = [bp for _, _, bp, _ in items if bp > 0]
    if not bps_nonzero:
        center = approx_size
    else:
        med = statistics.median(bps_nonzero)
        # trust approx_size if median looks off by >10%
        center = approx_size if abs(med - approx_size) / max(approx_size, 1) > 0.10 else int(med)

    def filter_by_window(pct):
        lo = int(center * (1 - pct))
        hi = int(center * (1 + pct))
        kept = [it for it in items if lo <= it[2] <= hi]
        return kept, lo, hi

    kept = items
    used_pct = None
    lo = hi = None
    for pct in pct_list:
        kept2, lo, hi = filter_by_window(pct)
        if len(kept2) >= args.min_keep:
            kept = kept2
            used_pct = pct
            break
    if used_pct is None:
        kept = items
        used_pct = -1
        lo = hi = None

    def rank_key(it):
        tig, fp, bp, depth = it
        return (abs(bp - approx_size) if bp > 0 else 10**12, -depth, -bp)

    kept_sorted = sorted(kept, key=rank_key)
    kept_final = kept_sorted[:max_keep]
    keep_names = set(tig for tig, *_ in kept_final)

    report_rows = []
    for tig, fp, bp, depth in sorted(items, key=lambda x: x[0]):
        keep = tig in keep_names
        if not keep:
            if used_pct != -1 and (bp < lo or bp > hi):
                reason = f"REJECT:length_outlier (window={lo}-{hi})"
            else:
                reason = "REJECT:ranked_out (cap_26)"
            fp.rename(reject_dir / fp.name)
        else:
            reason = "KEEP"
        report_rows.append((tig, bp, depth, "KEEP" if keep else "REJECT", reason))

    out = cluster_dir / "filter_report.tsv"
    with out.open("w") as f:
        f.write("tig\tbp\tdepth_x\tdecision\treason\n")
        for tig, bp, depth, decision, reason in report_rows:
            f.write(f"{tig}\t{bp}\t{depth:.3f}\t{decision}\t{reason}\n")

    print(f"[filter] approx_size={approx_size} center={center} pct={used_pct} kept={len(keep_names)} total={len(items)}")
    if missing:
        print("[filter] WARNING missing in contigs.phylip:", ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))

if __name__ == "__main__":
    main()
