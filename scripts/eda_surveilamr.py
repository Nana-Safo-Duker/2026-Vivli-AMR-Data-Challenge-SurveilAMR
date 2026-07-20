"""SurveilAMR exploratory data analysis — quick-look sample diagnostics.

Loads a fast in-memory sample of ATLAS (first 200k rows) for exploratory
plots that don't require the full chunked pipeline. Useful for rapid
iteration before running scripts/run_analysis.py on the complete 1M-row file.

Run: python scripts/eda_surveilamr.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "raw" / "atlas_vivli_2004_2024.csv"
OUT_DIR = ROOT / "outputs" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

INTERPRETATION_MAP = {
    "Susceptible": "S",
    "Resistant": "R",
    "Intermediate": "I",
    "S": "S",
    "I": "I",
    "R": "R",
}


def load_sample(nrows: int = 200_000) -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH, nrows=nrows, low_memory=False)
    interpret_cols = [c for c in df.columns if c.endswith("_I")]
    for col in interpret_cols:
        df[col] = df[col].map(INTERPRETATION_MAP)
    return df


def plot_top_species(df: pd.DataFrame) -> None:
    plt.figure(figsize=(10, 6))
    df["Species"].value_counts().head(10).plot(kind="bar", color="steelblue")
    plt.title("Top 10 Bacterial Species (Sample)")
    plt.ylabel("Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "eda_top_species.png")
    plt.close()


def plot_yearly_trend(df: pd.DataFrame, drug: str = "Meropenem") -> None:
    col = f"{drug}_I"
    if col not in df.columns:
        return
    ecoli = df[df["Species"] == "Escherichia coli"]
    trend = ecoli.groupby("Year")[col].value_counts(normalize=True).unstack()
    if trend.empty:
        return
    trend.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="RdYlGn_r")
    plt.title(f"E. coli {drug} Susceptibility Trend by Year (Sample)")
    plt.ylabel("Proportion")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"eda_ecoli_{drug.lower()}_trend.png")
    plt.close()


def plot_resistance_heatmap(df: pd.DataFrame) -> None:
    interpret_cols = [c for c in df.columns if c.endswith("_I")]
    top_species = df["Species"].value_counts().head(8).index
    subset = df[df["Species"].isin(top_species)]
    matrix = subset.groupby("Species")[interpret_cols].apply(lambda x: (x == "R").mean())
    matrix.columns = [c.replace("_I", "") for c in matrix.columns]
    plt.figure(figsize=(16, 6))
    sns.heatmap(matrix, cmap="Reds", annot=False, linewidths=0.3)
    plt.title("Resistance Heatmap: Top 8 Species (Sample)")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "eda_resistance_heatmap.png")
    plt.close()


def plot_country_resistance(df: pd.DataFrame) -> None:
    col = "Meropenem_I"
    kp = df[df["Species"] == "Klebsiella pneumoniae"]
    if col not in kp.columns or kp.empty:
        return
    country_trend = (
        kp.groupby("Country")[col]
        .apply(lambda x: (x == "R").mean())
        .sort_values(ascending=False)
        .head(15)
    )
    plt.figure(figsize=(10, 6))
    country_trend.plot(kind="barh", color="crimson")
    plt.title("K. pneumoniae Meropenem Resistance by Country (Sample, Top 15)")
    plt.xlabel("Resistance Proportion")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "eda_kp_meropenem_country.png")
    plt.close()


def main() -> None:
    sns.set(style="whitegrid")
    print(f"Loading sample from {DATA_PATH}...")
    df = load_sample()
    print(f"Sample size: {len(df):,} isolates")
    print(f"Columns: {len(df.columns)}")
    print(f"Years: {df['Year'].min()}\u2013{df['Year'].max()}")
    print(f"Countries: {df['Country'].nunique()}")

    plot_top_species(df)
    plot_yearly_trend(df, "Meropenem")
    plot_yearly_trend(df, "Ceftriaxone")
    plot_resistance_heatmap(df)
    plot_country_resistance(df)
    print(f"EDA figures saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
