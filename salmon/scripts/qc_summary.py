"""Resumo de QC por amostra a partir das saídas do Salmon.

Uso: python3 qc_summary.py samples.tsv <dir_salmon> qc_summary.tsv
"""
import csv
import json
import os
import sys

samples_tsv, salmon_dir, out = sys.argv[1:4]
with open(samples_tsv) as f:
    samples = list(csv.DictReader(f, delimiter="\t"))

cols = ["run", "condition", "mouse_id", "read_pairs_ena", "num_processed", "num_mapped",
        "percent_mapped", "library_type", "compatible_fragment_ratio", "salmon_version"]
with open(out, "w") as fout:
    fout.write("\t".join(cols) + "\n")
    for s in samples:
        d = os.path.join(salmon_dir, s["run"])
        meta = json.load(open(os.path.join(d, "aux_info", "meta_info.json")))
        lib = json.load(open(os.path.join(d, "lib_format_counts.json")))
        if meta["num_processed"] != int(s["read_pairs"]):
            print(f"AVISO {s['run']}: Salmon processou {meta['num_processed']} pares, "
                  f"ENA informa {s['read_pairs']}")
        row = [s["run"], s["condition"], s["mouse_id"], s["read_pairs"],
               meta["num_processed"], meta["num_mapped"], f"{meta['percent_mapped']:.2f}",
               lib["expected_format"], f"{lib['compatible_fragment_ratio']:.4f}",
               meta["salmon_version"]]
        fout.write("\t".join(map(str, row)) + "\n")

versions = set()
for s in samples:
    meta = json.load(open(os.path.join(salmon_dir, s["run"], "aux_info", "meta_info.json")))
    versions.add(meta["salmon_version"])
if len(versions) != 1:
    sys.exit(f"ERRO: amostras quantificadas com versões diferentes do Salmon: {versions}")
print(f"{len(samples)} amostras resumidas em {out}")
