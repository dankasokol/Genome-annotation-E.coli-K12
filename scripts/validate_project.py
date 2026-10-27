"""Read-only checks for the saved MG1655 annotation and its release.

Run from any directory: python /path/to/genome_annotation/scripts/validate_project.py
The script checks saved artifacts; it does not rerun annotation web services.
"""
from pathlib import Path
import hashlib
import json

from Bio import SeqIO
import nbformat
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(relative):
    return pd.read_csv(ROOT / relative, sep="\t")


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    genome = SeqIO.read(ROOT / "data/raw/genome.fna", "fasta")
    bakta = SeqIO.read(ROOT / "results/structural/bakta/ecoli.gbff", "genbank")
    ncbi_path = ROOT / "data/processed/NC_000913.3.gbk"
    require(str(genome.seq).upper() == str(bakta.seq).upper(), "Bakta sequence differs from input")
    require(len(genome.seq) == 4_641_652, "Unexpected input genome length")
    cds = [f for f in bakta.features if f.type == "CDS"]
    ids = [f.qualifiers["locus_tag"][0] for f in cds]
    require(len(cds) == len(set(ids)) == 4315, "CDS count or locus tags changed")
    protein_records = list(SeqIO.parse(ROOT / "results/structural/bakta/ecoli.faa", "fasta"))
    proteins = {p.id: str(p.seq).upper().rstrip("*") for p in protein_records}
    require(len(protein_records) == len(proteins) == len(cds), "Duplicated or missing protein IDs")
    require(set(proteins) == set(ids), "FASTA and GenBank IDs differ")
    for f in cds:
        key = f.qualifiers["locus_tag"][0]
        require(proteins[key] == f.qualifiers["translation"][0].upper().rstrip("*"), f"Translation differs: {key}")

    integrated = table("results/integrated/ecoli_complete_annotation.tsv")
    primary = table("results/release_v1_bakta/tables/bakta_CDS_functional_annotation.tsv")
    for df in (integrated, primary):
        require(len(df) == 4315 and df["query"].is_unique, "Invalid functional catalog grain")
        require(set(df["query"]) == set(ids), "Functional catalog and Bakta IDs differ")
        lookup = df.set_index("query")
        for f in cds:
            r = lookup.loc[f.qualifiers["locus_tag"][0]]
            require((int(r.start), int(r.end), int(r.strand)) ==
                    (int(f.location.start) + 1, int(f.location.end), f.location.strand),
                    f"Coordinates differ: {r.name}")
    pd.testing.assert_frame_equal(integrated.sort_values("query").reset_index(drop=True),
                                  primary[integrated.columns].sort_values("query").reset_index(drop=True))
    # Keep the #query header, while excluding the ## metadata lines.
    from io import StringIO
    lines = (ROOT / "results/functional/eggnog/query.emapper.annotations").read_text().splitlines()
    eggnog = pd.read_csv(StringIO("\n".join(l for l in lines if not l.startswith("##"))), sep="\t", na_values=["-"])
    eggnog.columns = eggnog.columns.str.lstrip("#")
    require(eggnog["query"].is_unique and set(eggnog["query"]) <= set(ids), "Unexpected eggNOG IDs")
    pd.testing.assert_frame_equal(
        integrated.set_index("query").loc[eggnog["query"], eggnog.columns.drop("query")].reset_index(),
        eggnog.reset_index(drop=True), check_dtype=False)

    hypothetical = integrated[integrated.bakta_product.eq("hypothetical protein")]
    remaining = hypothetical[hypothetical[["GOs", "KEGG_ko", "PFAMs"]].isna().all(axis=1)]
    selected = list(SeqIO.parse(ROOT / "results/functional/interpro/hypothetical.faa", "fasta"))
    require({r.id for r in selected} == set(hypothetical["query"]), "InterPro input set differs")
    md5_to_ids = {}
    for p in selected:
        md5_to_ids.setdefault(hashlib.md5(str(p.seq).upper().encode("ascii")).hexdigest().upper(), []).append(p.id)
    batch_records = [r for f in sorted((ROOT / "results/functional/interpro").glob("hypothetical_batch_*.faa"))
                     for r in SeqIO.parse(f, "fasta")]
    require([(r.id, str(r.seq)) for r in selected] == [(r.id, str(r.seq)) for r in batch_records], "InterPro batches differ")
    raw = pd.concat([pd.read_csv(f, sep="\t", header=None, dtype=str, na_values=["-"])
                     for f in sorted((ROOT / "results/functional/interpro").glob("result_interPro_*.tsv"))], ignore_index=True)
    require(set(raw[1].str.upper()) <= set(md5_to_ids), "Unknown InterPro sequence digest")
    require((raw[2].astype(int) == raw[1].str.upper().map(
        {h: len(proteins[queries[0]]) for h, queries in md5_to_ids.items()})).all(), "InterPro lengths differ")
    require(not integrated.loc[~integrated["query"].isin(hypothetical["query"]), "signatures"].notna().any(),
            "InterPro results assigned to unsubmitted proteins")
    integrated_lookup = integrated.set_index("query")
    for digest, queries in md5_to_ids.items():
        rows = raw[raw[1].str.upper().eq(digest)]
        for source_col, result_col in [(4, "signatures"), (11, "interpro_ids"), (12, "interpro_domains"), (13, "interpro_go"), (3, "analyses")]:
            expected = "; ".join(sorted({str(v) for v in rows[source_col].dropna() if str(v).strip()}))
            for query in queries:
                actual = integrated_lookup.at[query, result_col]
                require(("" if pd.isna(actual) else actual) == expected, f"InterPro aggregation differs: {query}, {result_col}")

    mapping = table("results/release_v1_bakta/tables/defense_genes_bakta_mapping.tsv")
    systems = table("results/release_v1_bakta/tables/defense_systems.tsv")
    phastest = table("results/release_v1_bakta/tables/phastest_regions.tsv")
    defense_proteins = {r.id: str(r.seq).upper().rstrip("*") for r in SeqIO.parse(
        ROOT / "results/specialized/defensefinder/defensefinder_component_proteins.faa", "fasta")}
    require(len(mapping) == 23 and len(systems) == 9 and len(phastest) == 5, "Specialized counts changed")
    require(set(mapping.sys_id) <= set(systems.sys_id), "Unknown defense system")
    lookup = integrated.set_index("query")
    for r in mapping.itertuples():
        if r.mapping_status == "exact_coordinates_and_protein":
            b = lookup.loc[r.bakta_query]
            require((r.start_1based, r.end_1based, {"+": 1, "-": -1}[r.strand]) ==
                    (b.start, b.end, b.strand), "DefenseFinder coordinate mismatch")
            require(defense_proteins[r.hit_id] == proteins[r.bakta_query], "DefenseFinder protein mismatch")
    for df in (mapping, systems, phastest):
        require(((df.start_1based >= 1) & (df.end_1based <= len(genome)) &
                 (df.end_1based >= df.start_1based)).all(), "Invalid specialized interval")
    for r in systems.itertuples():
        components = mapping[mapping.sys_id.eq(r.sys_id)]
        require((components.start_1based.min(), components.end_1based.max()) ==
                (r.start_1based, r.end_1based), "System interval differs from components")
    require(((phastest.end_1based - phastest.start_1based + 1) == phastest.length_bp).all(), "PHASTEST length mismatch")

    gff = ROOT / "results/structural/bakta/ecoli.gff3"
    seqids = set()
    with gff.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip().split("\t")
                require(len(fields) == 9, "Invalid Bakta GFF row")
                seqids.add(fields[0])
    require(len(seqids) == 1, "Expected one GFF sequence")
    for name, df in [("phastest_regions", phastest), ("defense_systems", systems)]:
        path = ROOT / f"results/release_v1_bakta/tracks/{name}.gff3"
        rows = [l.split("\t") for l in path.read_text().splitlines() if l and not l.startswith("#")]
        require(len(rows) == len(df), "GFF track count mismatch")
        require({r[0] for r in rows} == seqids and all(len(r) == 9 for r in rows), "Invalid track sequence or format")
        require([(int(r[3]), int(r[4])) for r in rows] == list(zip(df.start_1based, df.end_1based)), "Track intervals differ")

    summary = json.loads((ROOT / "results/release_v1_bakta/qc/summary.json").read_text())
    require(summary["primary_bakta_cds"] == len(cds), "Summary CDS count differs")
    require(summary["primary_bakta_features_total"] == len(bakta.features), "Summary feature count differs")
    require(summary["hypothetical_proteins_Bakta"] == len(hypothetical), "Summary hypothetical count differs")
    require(summary["DefenseFinder_exact_matches_to_Bakta"] == int(mapping.mapping_status.eq("exact_coordinates_and_protein").sum()), "Summary defense count differs")
    if ncbi_path.is_file():
        ncbi = SeqIO.read(ncbi_path, "genbank")
        require(str(genome.seq).upper() == str(ncbi.seq).upper(), "NCBI sequence differs")
        ref_cds = [f for f in ncbi.features if f.type == "CDS"]
        nc_locations = {str(f.location) for f in ref_cds}
        exact = sum(str(f.location) in nc_locations for f in cds)
        require(summary["NCBI_validation_only"]["exact_cds_coordinate_matches"] == exact, "Summary NCBI matches differ")
    checks = table("results/release_v1_bakta/qc/checksums.tsv")
    for row in checks.itertuples():
        require(sha256(ROOT / row.path) == row.sha256, f"Checksum differs: {row.path}")
    for suffix in ("gbff", "gff3"):
        require(sha256(ROOT / f"results/structural/bakta/ecoli.{suffix}") ==
                sha256(ROOT / f"results/release_v1_bakta/ecoli_bakta_primary.{suffix}"), "Release copy changed")

    notebooks = sorted((ROOT / "notebooks").glob("*.ipynb"))
    for path in notebooks:
        nb = nbformat.read(path, as_version=4)
        nbformat.validate(nb)
        for index, cell in enumerate(nb.cells):
            if cell.cell_type == "code":
                compile(cell.source, f"{path.name}:{index}", "exec")
                require(cell.execution_count is not None, f"Unexecuted cell: {path.name}:{index}")
                require(not any(o.output_type == "error" for o in cell.outputs), f"Saved error: {path.name}:{index}")
    print(json.dumps({"status": "PASS", "genome_length_bp": len(genome), "bakta_cds": len(cds),
                      "eggnog_rows": len(eggnog), "hypothetical_proteins": len(hypothetical),
                      "without_eggnog_GO_KO_Pfam": len(remaining),
                      "remaining_with_IPR": int(remaining.interpro_ids.notna().sum()),
                      "remaining_with_any_signature": int(remaining.signatures.notna().sum()),
                      "defense_components": len(mapping), "phastest_regions": len(phastest),
                      "checksums_verified": len(checks), "notebooks_checked": len(notebooks)}, indent=2))


if __name__ == "__main__":
    main()
