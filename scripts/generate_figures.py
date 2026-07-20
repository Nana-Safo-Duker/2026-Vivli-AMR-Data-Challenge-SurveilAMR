"""Generate the primary ATLAS publication figures for SurveilAMR.

Reads data/processed/dataset_summary.json and the yearly/country resistance
CSVs produced by scripts/run_analysis.py, and writes 7 figures to
outputs/figures/ (fig1 .. fig6, fig10).

Run: python scripts/generate_figures.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
PROC_DIR = ROOT / "data" / "processed"
OUT_DIR = ROOT / "outputs" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"figure.dpi": 150, "font.size": 10})
sns.set_theme(style="whitegrid", palette="muted")


def load_summary() -> dict:
    with open(PROC_DIR / "dataset_summary.json", encoding="utf-8") as f:
        return json.load(f)


def fig_species_distribution(summary: dict) -> None:
    species = summary["top_species"]
    df = pd.DataFrame(list(species.items()), columns=["Species", "Count"])
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=df, x="Count", y="Species", hue="Species", legend=False, ax=ax, palette="Blues_r")
    ax.set_title("Top 15 Pathogens in ATLAS Antibiotics Surveillance (2004\u20132024)")
    ax.set_xlabel("Number of Isolates")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_species_distribution.png", bbox_inches="tight")
    plt.close(fig)


def fig_country_distribution(summary: dict) -> None:
    countries = summary["top_countries"]
    df = pd.DataFrame(list(countries.items()), columns=["Country", "Count"])
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=df, x="Count", y="Country", hue="Country", legend=False, ax=ax, palette="Greens_r")
    ax.set_title("Top 15 Contributing Countries (ATLAS)")
    ax.set_xlabel("Number of Isolates")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig2_country_distribution.png", bbox_inches="tight")
    plt.close(fig)


def fig_resistance_trends() -> None:
    df = pd.read_csv(PROC_DIR / "surveilamr_resistance_by_year.csv")
    focus = [
        ("Escherichia coli", "Meropenem"),
        ("Klebsiella pneumoniae", "Meropenem"),
        ("Pseudomonas aeruginosa", "Meropenem"),
        ("Staphylococcus aureus", "Oxacillin"),
        ("Acinetobacter baumannii", "Meropenem"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
    axes = axes.flatten()
    for i, (sp, ab) in enumerate(focus):
        sub = df[(df["Species"] == sp) & (df["Antibiotic"] == ab)].sort_values("Year")
        if sub.empty:
            axes[i].set_visible(False)
            continue
        axes[i].plot(sub["Year"], sub["Resistance_Pct"], marker="o", color="#c0392b", linewidth=2)
        axes[i].set_title(f"{sp.split()[0]} \u2013 {ab}", fontsize=9)
        axes[i].set_xlabel("Year")
        axes[i].set_ylabel("Resistance (%)")
        axes[i].set_ylim(0, max(100, sub["Resistance_Pct"].max() * 1.1))
    axes[-1].set_visible(False)
    fig.suptitle("Temporal Trends in Key Antibiotic Resistance (ESKAPE Pathogens)", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_resistance_trends.png", bbox_inches="tight")
    plt.close(fig)


def fig_heatmap_country_resistance() -> None:
    df = pd.read_csv(PROC_DIR / "surveilamr_resistance_by_country.csv")
    top_countries = df.groupby("Country")["Total_Isolates"].sum().nlargest(12).index
    sub = df[df["Country"].isin(top_countries)]
    pivot = sub.pivot_table(index="Species", columns="Country", values="Resistance_Pct", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(14, 5))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd", linewidths=0.5, ax=ax)
    ax.set_title("Country-Level Resistance to Primary Antibiotics (Priority Pathogens)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_country_resistance_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def fig_gene_prevalence(summary: dict) -> None:
    genes = summary["gene_detection_pct"]
    df = pd.DataFrame(list(genes.items()), columns=["Gene", "Detection_Pct"])
    df = df.sort_values("Detection_Pct", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=df, x="Detection_Pct", y="Gene", hue="Gene", legend=False, ax=ax, palette="Purples_r")
    ax.set_title("Beta-Lactamase Gene Detection in ATLAS Isolates")
    ax.set_xlabel("Detection Rate (%)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_gene_prevalence.png", bbox_inches="tight")
    plt.close(fig)


def fig_mdr_proxy(summary: dict) -> None:
    mdr = summary["multidrug_resistance_proxy"]
    df = pd.DataFrame(list(mdr.items()), columns=["Species", "MDR_Proxy_Pct"])
    df["Label"] = df["Species"].str.split().str[0]
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.barplot(data=df, x="Label", y="MDR_Proxy_Pct", hue="Label", legend=False, ax=ax, palette="Oranges")
    ax.set_title("Multi-Drug Resistance Proxy (\u2265 2 Resistant Classes)")
    ax.set_ylabel("Proportion (%)")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig6_mdr_proxy.png", bbox_inches="tight")
    plt.close(fig)


def fig_africa_comparison(summary: dict) -> None:
    rows = []
    for sp, entry in summary["resistance_trends"].items():
        years = entry.get("yearly_resistance_pct", {})
        if not years:
            continue
        latest_year = max(int(y) for y in years)
        rows.append(
            {
                "Species": sp.split()[0],
                "Global_latest": years[str(latest_year)],
                "Africa": entry.get("africa_resistance_pct"),
            }
        )
    df = pd.DataFrame(rows).dropna()
    if df.empty:
        return
    melted = df.melt(id_vars="Species", var_name="Region", value_name="Resistance_Pct")
    fig, ax = plt.subplots(figsize=(9, 4.5))
    sns.barplot(data=melted, x="Species", y="Resistance_Pct", hue="Region", ax=ax)
    ax.set_title("Primary-Antibiotic Resistance: Global Latest Year vs Africa (ATLAS)")
    ax.set_ylabel("Resistance (%)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig10_africa_vs_global.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    summary = load_summary()
    fig_species_distribution(summary)
    fig_country_distribution(summary)
    fig_resistance_trends()
    fig_heatmap_country_resistance()
    fig_gene_prevalence(summary)
    fig_mdr_proxy(summary)
    fig_africa_comparison(summary)
    print(f"Figures saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
