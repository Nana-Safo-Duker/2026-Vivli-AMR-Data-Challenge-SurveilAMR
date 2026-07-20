"""SurveilAMR — KEYSTONE, DREAM, and GASAR (Study III) analysis.

Produces per-dataset summaries (data/processed/supplementary_datasets_summary.json),
tidy CSV extracts, and matplotlib/seaborn figures for the omadacycline (KEYSTONE),
bedaquiline (DREAM), and Gram-negative mechanism (GASAR) surveillance datasets.

Run: python scripts/analyze_supplementary.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "outputs" / "figures"
PROC.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

DREAM_PATH = RAW / "BEDAQUILINE DREAM DATASET FOR VIVLI - 06-06-2022.xlsx"
GASAR_PATH = RAW / "GASAR Study III (n=494)_updated.xlsx"
KEYSTONE_PATH = RAW / "Omadacycline_2015 to 2025_Surveillance_data.xlsx"

plt.rcParams.update({"figure.dpi": 150, "font.size": 10})
sns.set_theme(style="whitegrid")

MIC_STRIP_RE = re.compile(r"[\u2264\u2265<>]")  # strips ≤ ≥ < >


def parse_mic(value) -> float:
    """Convert MIC strings such as '<=0.03' or '>32' to a numeric value."""
    if pd.isna(value):
        return np.nan
    text = MIC_STRIP_RE.sub("", str(value)).strip()
    try:
        return float(text)
    except ValueError:
        return np.nan


def normalize_tb_subtype(value) -> str:
    """Collapse the ~40 free-text DREAM SubType spellings into 4 clean buckets."""
    s = str(value).upper().replace("\xa0", " ").strip()
    s = re.sub(r"\s+", " ", s)
    s_nospace = s.replace(" ", "").replace("-", "")
    if "PRE" in s_nospace and "XDR" in s_nospace:
        return "Pre-XDR"
    if "XDR" in s_nospace:
        return "XDR"
    if "MDR" in s_nospace:
        return "MDR"
    if "RR" in s_nospace or "RIF" in s_nospace or "RMP" in s_nospace or s in {"R", "*R"}:
        return "RR/RIF-resistant"
    return "Other/Unclassified"


def analyze_dream() -> dict:
    df = pd.read_excel(DREAM_PATH)
    df["Continent"] = df["Continent"].astype(str).str.strip().str.title()
    df["Country"] = df["Country"].astype(str).str.strip()
    df["BDQ_MIC"] = df["BDQ Broth"].map(parse_mic)
    df["Subtype_Norm"] = df["SubType"].map(normalize_tb_subtype)

    summary = {
        "n": int(len(df)),
        "years": [int(df["Year Collected"].min()), int(df["Year Collected"].max())],
        "n_countries": int(df["Country"].nunique()),
        "continents": df["Continent"].value_counts().to_dict(),
        "countries": df["Country"].value_counts().to_dict(),
        "subtype_normalized": df["Subtype_Norm"].value_counts().to_dict(),
        "bdq_mic_median": float(df["BDQ_MIC"].median()),
        "bdq_mic_mean": round(float(df["BDQ_MIC"].mean()), 4),
        "bdq_gt_0_12_pct": round(float((df["BDQ_MIC"] > 0.12).mean() * 100), 2),
        "bdq_gt_0_25_pct": round(float((df["BDQ_MIC"] > 0.25).mean() * 100), 2),
        "africa_n": int((df["Continent"] == "Africa").sum()),
        "bdq_median_by_continent": df.groupby("Continent")["BDQ_MIC"].median().round(3).to_dict(),
        "bdq_median_by_year": {int(k): v for k, v in df.groupby("Year Collected")["BDQ_MIC"].median().round(3).to_dict().items()},
        "bdq_median_by_subtype": df.groupby("Subtype_Norm")["BDQ_MIC"].median().round(3).to_dict(),
    }

    yearly = (
        df.groupby("Year Collected")
        .agg(n=("BDQ_MIC", "size"), median_bdq=("BDQ_MIC", "median"), pct_gt_0_12=("BDQ_MIC", lambda x: (x > 0.12).mean() * 100))
        .reset_index()
    )
    yearly.to_csv(PROC / "dream_bdq_by_year.csv", index=False)

    cont = (
        df.groupby("Continent")
        .agg(n=("BDQ_MIC", "size"), median_bdq=("BDQ_MIC", "median"), mean_bdq=("BDQ_MIC", "mean"))
        .reset_index()
    )
    cont.to_csv(PROC / "dream_bdq_by_continent.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    order = cont.sort_values("n", ascending=False)["Continent"]
    sns.barplot(data=cont, x="Continent", y="median_bdq", order=order, ax=axes[0], color="#2980b9")
    axes[0].set_title("DREAM: Median Bedaquiline MIC by Continent")
    axes[0].set_ylabel("Median BDQ broth MIC (\u00b5g/mL)")
    axes[0].tick_params(axis="x", rotation=20)

    sns.lineplot(data=yearly, x="Year Collected", y="median_bdq", marker="o", ax=axes[1], color="#c0392b")
    axes[1].set_title("DREAM: Median BDQ MIC Over Time")
    axes[1].set_ylabel("Median BDQ broth MIC (\u00b5g/mL)")
    fig.tight_layout()
    fig.savefig(FIG / "fig7_dream_bdq.png", bbox_inches="tight")
    plt.close(fig)

    return summary


def analyze_gasar() -> dict:
    df = pd.read_excel(GASAR_PATH, sheet_name="Sheet1")
    df["Country"] = df["Country"].astype(str).str.strip()
    df["Poly_MIC"] = pd.to_numeric(df["Polymyxin B MIC (mcg/ml)"], errors="coerce")

    def phenotype_bucket(val) -> str:
        s = str(val).upper()
        if "NON ESBL" in s and "NON MBL" in s:
            return "Non-ESBL/Non-MBL"
        has_esbl = "ESBL" in s
        has_mbl = "MBL" in s
        has_carb = "CARBAPENEMASE" in s
        if has_esbl and has_mbl:
            return "ESBL+MBL"
        if has_mbl:
            return "MBL"
        if has_esbl:
            return "ESBL"
        if has_carb:
            return "Carbapenemase"
        return "Other"

    df["Phenotype_Bucket"] = df["Phenotypic Combination"].map(phenotype_bucket)

    summary = {
        "n": int(len(df)),
        "years": [int(df["Year"].min()), int(df["Year"].max())],
        "countries": df["Country"].value_counts().to_dict(),
        "species": df["Species"].value_counts().to_dict(),
        "phenotype_buckets": df["Phenotype_Bucket"].value_counts().to_dict(),
        "gene_top": df["Gene Combination"].value_counts().head(12).to_dict(),
        "poly_mic_median": float(df["Poly_MIC"].median()),
        "poly_mic_ge4_pct": round(float((df["Poly_MIC"] >= 4).mean() * 100), 2),
        "mbl_or_carb_pct": round(
            float(df["Phenotype_Bucket"].isin(["MBL", "ESBL+MBL", "Carbapenemase"]).mean() * 100), 2
        ),
    }

    pheno = df["Phenotype_Bucket"].value_counts().rename_axis("Phenotype").reset_index(name="Count")
    pheno.to_csv(PROC / "gasar_phenotype_counts.csv", index=False)

    genes = df["Gene Combination"].value_counts().head(12).rename_axis("Gene").reset_index(name="Count")
    genes.to_csv(PROC / "gasar_gene_counts.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.barplot(data=pheno, x="Count", y="Phenotype", hue="Phenotype", legend=False, ax=axes[0], palette="Purples_r")
    axes[0].set_title("GASAR: Phenotypic Resistance Classes (India)")
    genes_plot = genes[genes["Gene"].astype(str) != "-"]
    sns.barplot(data=genes_plot, x="Count", y="Gene", hue="Gene", legend=False, ax=axes[1], palette="Greens_r")
    axes[1].set_title("GASAR: Top Detected Gene Combinations")
    fig.tight_layout()
    fig.savefig(FIG / "fig8_gasar_mechanisms.png", bbox_inches="tight")
    plt.close(fig)

    return summary


def analyze_keystone() -> dict:
    df = pd.read_excel(KEYSTONE_PATH)
    df.columns = [str(c).replace("\n", " ").strip() for c in df.columns]
    df["Oma_MIC"] = df["Omadacycline"].map(parse_mic)

    priority = [
        "Staphylococcus aureus",
        "Escherichia coli",
        "Klebsiella pneumoniae",
        "Streptococcus pneumoniae",
        "Enterococcus faecalis",
        "Pseudomonas aeruginosa",
    ]

    mic_rows = []
    for sp in priority:
        sub = df.loc[df["Organism"] == sp, "Oma_MIC"].dropna()
        if len(sub) == 0:
            continue
        mic_rows.append(
            {
                "Organism": sp,
                "n": int(len(sub)),
                "MIC50": float(sub.quantile(0.5)),
                "MIC90": float(sub.quantile(0.9)),
                "mean": round(float(sub.mean()), 4),
            }
        )
    mic_df = pd.DataFrame(mic_rows)
    mic_df.to_csv(PROC / "keystone_omadacycline_mic.csv", index=False)

    yearly = (
        df[df["Organism"].isin(priority)]
        .groupby(["Study Year", "Organism"])["Oma_MIC"]
        .median()
        .reset_index(name="MIC50")
    )
    yearly.to_csv(PROC / "keystone_omadacycline_by_year.csv", index=False)

    africa = df[df["Continent"].astype(str).str.contains("Africa", case=False, na=False)]
    summary = {
        "n": int(len(df)),
        "years": [int(df["Study Year"].min()), int(df["Study Year"].max())],
        "n_species": int(df["Organism"].nunique()),
        "n_countries": int(df["Country"].nunique()),
        "top_species": df["Organism"].value_counts().head(12).to_dict(),
        "top_countries": df["Country"].value_counts().head(12).to_dict(),
        "continents": df["Continent"].value_counts().to_dict(),
        "africa_n": int(len(africa)),
        "africa_countries": africa["Country"].value_counts().to_dict() if len(africa) else {},
        "omadacycline_mic_overall_median": float(df["Oma_MIC"].median()),
        "mic_summary": mic_rows,
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.barplot(data=mic_df, x="Organism", y="MIC90", hue="Organism", legend=False, ax=axes[0], palette="Oranges_r")
    axes[0].set_title("KEYSTONE: Omadacycline MIC90 by Priority Pathogen")
    axes[0].set_ylabel("MIC90 (\u00b5g/mL)")
    axes[0].tick_params(axis="x", rotation=30)
    for label in axes[0].get_xticklabels():
        label.set_ha("right")

    focus = yearly[yearly["Organism"].isin(["Staphylococcus aureus", "Escherichia coli", "Klebsiella pneumoniae"])]
    sns.lineplot(data=focus, x="Study Year", y="MIC50", hue="Organism", marker="o", ax=axes[1])
    axes[1].set_title("KEYSTONE: Omadacycline MIC50 Trends")
    axes[1].set_ylabel("MIC50 (\u00b5g/mL)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig9_keystone_omadacycline.png", bbox_inches="tight")
    plt.close(fig)

    return summary


def main() -> None:
    summary = {
        "dream": analyze_dream(),
        "gasar": analyze_gasar(),
        "keystone": analyze_keystone(),
    }
    with open(PROC / "supplementary_datasets_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)
    print("KEYSTONE n =", summary["keystone"]["n"], "| DREAM n =", summary["dream"]["n"], "| GASAR n =", summary["gasar"]["n"])
    print(f"Saved summaries and figures to {PROC} and {FIG}")


if __name__ == "__main__":
    main()
