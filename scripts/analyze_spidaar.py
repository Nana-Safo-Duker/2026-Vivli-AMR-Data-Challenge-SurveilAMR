"""SurveilAMR — SPIDAAR RWE Study analysis (Pfizer-funded real-world evidence study).

SPIDAAR ("Surveillance of Pathogens and Drug resistance In an African Antimicrobial
Resistance network") contributes hospital-level isolate (n=244) and patient (n=336)
data across 4 sub-Saharan African countries (Ghana, Kenya, Malawi, Uganda).

This script decodes the numeric survey codes using the accompanying codebook
sheets, then produces tidy CSV extracts and figures describing:
  * organism-group distribution by country
  * multidrug resistance (MDR), MRSA, and third-generation cephalosporin
    resistance (3GC-R) flag rates among tested isolates
  * patient demographics, healthcare-associated infection (HAI) categories,
    device exposure, and in-hospital mortality by disease severity

Run: python scripts/analyze_spidaar.py
"""
from __future__ import annotations

import json
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

ISOLATE_PATH = RAW / "spidaar_isolatedata.xls"
PATIENT_PATH = RAW / "spidaar_patientdata.xls"

plt.rcParams.update({"figure.dpi": 150, "font.size": 10})
sns.set_theme(style="whitegrid")

GROUP_LABELS = {
    1: "Staphylococcus aureus", 2: "Streptococcus", 3: "Enterobacterales",
    4: "Pseudomonads", 5: "Acinetobacters", 6: "Enterococcus",
    7: "Other Staphylococci", 8: "Other Gram-negative", 9: "Other Gram-positive",
}

CHAICAT_ISOLATE_LABELS = {
    1: "BSI", 2: "cUTI", 3: "cSSTI", 4: "HAP", 5: "cIAI", 6: "BSI + cUTI",
    7: "BSI + cSSTI", 8: "HAP + cSSTI", 9: "BSI + cUTI + cIAI", 10: "No HAI confirmed",
}

WARD_LABELS = {1: "ICU", 2: "Internal medicine", 3: "Surgery wards", 4: "Other wards"}

AGEGR_LABELS = {
    0: "<1", 1: "1-5", 2: "6-10", 3: "11-15", 4: "16-20", 5: "21-25", 6: "26-30",
    7: "31-35", 8: "36-40", 9: "41-45", 10: "46-50", 11: "51-55", 12: "56-60",
    13: "61-65", 14: "65+",
}

SEX_LABELS = {0: "Male", 1: "Female", 3: "Missing"}
DISEV_LABELS = {1: "Mild", 2: "Moderate", 3: "Severe"}
DEAD_LABELS = {0: "Alive", 1: "Deceased", 9: "Deceased (censored)"}


def analyze_isolates() -> dict:
    df = pd.read_excel(ISOLATE_PATH, sheet_name="data")
    df["Group_Label"] = df["group"].map(GROUP_LABELS).fillna("Unclassified")
    df["ctry"] = df["ctry"].astype(str).str.strip()

    def flag_rate(col: str, unknown_codes: set[int]) -> dict:
        vals = df[col][~df[col].isin(unknown_codes)]
        n = len(vals)
        if n == 0:
            return {"n_tested": 0, "positive_pct": None}
        return {"n_tested": int(n), "positive_pct": round(float((vals == 1).mean() * 100), 2)}

    flags = {
        "MDR": flag_rate("mdr", {99}),
        "MRSA": flag_rate("mrsa", {99}),
        "3GC_resistant": flag_rate("c3r", {9}),
    }

    by_country_group = (
        df.groupby(["ctry", "Group_Label"]).size().rename("Count").reset_index()
    )
    by_country_group.to_csv(PROC / "spidaar_isolate_groups_by_country.csv", index=False)

    flags_df = pd.DataFrame(
        [{"Flag": k, **v} for k, v in flags.items()]
    )
    flags_df.to_csv(PROC / "spidaar_amr_flags.csv", index=False)

    summary = {
        "n": int(len(df)),
        "countries": df["ctry"].value_counts().to_dict(),
        "group_distribution": df["Group_Label"].value_counts().to_dict(),
        "specimen_type_top": df["stype"].value_counts().head(10).to_dict(),
        "clinical_relevance": df["clinrel"].value_counts().to_dict(),
        "hai_category": {
            CHAICAT_ISOLATE_LABELS.get(k, str(k)): int(v)
            for k, v in df["chaicat"].value_counts().to_dict().items()
        },
        "amr_flags": flags,
    }

    top_groups = df["Group_Label"].value_counts().head(6).index
    plot_df = by_country_group[by_country_group["Group_Label"].isin(top_groups)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    pivot = plot_df.pivot_table(index="ctry", columns="Group_Label", values="Count", aggfunc="sum").fillna(0)
    pivot.plot(kind="bar", stacked=True, ax=axes[0], colormap="tab10")
    axes[0].set_title("SPIDAAR: Isolate Organism Groups by Country")
    axes[0].set_xlabel("Country")
    axes[0].set_ylabel("Isolates")
    axes[0].legend(fontsize=7, title=None)
    axes[0].tick_params(axis="x", rotation=0)

    sns.barplot(data=flags_df.dropna(subset=["positive_pct"]), x="Flag", y="positive_pct", hue="Flag", legend=False, ax=axes[1], palette="Reds_r")
    axes[1].set_title("SPIDAAR: Resistance Flag Positivity Among Tested Isolates")
    axes[1].set_ylabel("Positive (%)")
    for i, row in flags_df.dropna(subset=["positive_pct"]).reset_index(drop=True).iterrows():
        axes[1].text(i, row["positive_pct"] + 1.5, f"n={row['n_tested']}", ha="center", fontsize=8)
    axes[1].set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(FIG / "fig11_spidaar_isolate_groups.png", bbox_inches="tight")
    plt.close(fig)

    return summary


def analyze_patients() -> dict:
    df = pd.read_excel(PATIENT_PATH, sheet_name="data")
    df["ctry"] = df["ctry"].astype(str).str.strip()
    df["Ward_Label"] = df["ward"].map(WARD_LABELS).fillna("Unknown")
    df["Age_Label"] = df["agegr"].map(AGEGR_LABELS).fillna("Unknown")
    df["Sex_Label"] = df["sex"].map(SEX_LABELS).fillna("Unknown")
    df["Disease_Severity"] = df["disev"].map(DISEV_LABELS).fillna("Unknown")
    df["Outcome"] = df["dead"].map(DEAD_LABELS).fillna("Unknown")
    df["HAI_Category"] = df["chaicat"].map(CHAICAT_ISOLATE_LABELS).fillna("Unclassified")

    los = pd.to_numeric(df["los"], errors="coerce")

    demographics = {
        "n": int(len(df)),
        "countries": df["ctry"].value_counts().to_dict(),
        "sex": df["Sex_Label"].value_counts().to_dict(),
        "age_group": df["Age_Label"].value_counts().to_dict(),
        "ward": df["Ward_Label"].value_counts().to_dict(),
        "disease_severity": df["Disease_Severity"].value_counts().to_dict(),
        "outcome": df["Outcome"].value_counts().to_dict(),
        "hai_category": df["HAI_Category"].value_counts().to_dict(),
        "device_use_pct": round(float((df["devyn"] == 1).mean() * 100), 2),
        "surgery_pct": round(float((df["surgyn"] == 1).mean() * 100), 2),
        "los_days_median": float(los.median()),
        "los_days_mean": round(float(los.mean()), 2),
        "mortality_pct": round(float((df["dead"] == 1).mean() * 100), 2),
        "pathogen_positive_pct": round(float((df["pathp"] == 1).mean() * 100), 2),
    }

    los_by_outcome = df.assign(LOS=los).groupby("Outcome")["LOS"].agg(["count", "median", "mean"]).reset_index()
    los_by_outcome.to_csv(PROC / "spidaar_los_by_outcome.csv", index=False)

    severity_outcome = (
        df.groupby(["Disease_Severity", "Outcome"]).size().rename("Count").reset_index()
    )
    severity_outcome.to_csv(PROC / "spidaar_severity_by_outcome.csv", index=False)

    demo_df = pd.DataFrame(
        [{"Category": "Age group", "Value": k, "Count": v} for k, v in demographics["age_group"].items()]
        + [{"Category": "Ward", "Value": k, "Count": v} for k, v in demographics["ward"].items()]
    )
    demo_df.to_csv(PROC / "spidaar_patient_demographics.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ward_order = ["ICU", "Internal medicine", "Surgery wards", "Other wards", "Unknown"]
    ward_counts = df["Ward_Label"].value_counts().reindex(ward_order).dropna()
    sns.barplot(x=ward_counts.index, y=ward_counts.values, hue=ward_counts.index, legend=False, ax=axes[0], palette="Blues_r")
    axes[0].set_title("SPIDAAR: Patients by Hospital Ward")
    axes[0].set_ylabel("Patients")
    axes[0].tick_params(axis="x", rotation=20)

    sev_order = ["Mild", "Moderate", "Severe", "Unknown"]
    mortality_by_sev = (
        df[df["Outcome"].isin(["Alive", "Deceased"])]
        .assign(Disease_Severity=lambda d: d["Disease_Severity"])
        .groupby("Disease_Severity")["Outcome"]
        .apply(lambda x: (x == "Deceased").mean() * 100)
        .reindex(sev_order)
        .dropna()
    )
    sns.barplot(x=mortality_by_sev.index, y=mortality_by_sev.values, hue=mortality_by_sev.index, legend=False, ax=axes[1], palette="Reds")
    axes[1].set_title("SPIDAAR: In-Hospital Mortality by Disease Severity")
    axes[1].set_ylabel("Mortality (%)")

    los_plot = df.assign(LOS=los).dropna(subset=["LOS"])
    sns.boxplot(data=los_plot, x="Outcome", y="LOS", hue="Outcome", legend=False, ax=axes[2], palette="Set2")
    axes[2].set_title("SPIDAAR: Length of Stay by Outcome")
    axes[2].set_ylabel("Length of stay (days)")

    fig.tight_layout()
    fig.savefig(FIG / "fig12_spidaar_patient_outcomes.png", bbox_inches="tight")
    plt.close(fig)

    return demographics


def main() -> None:
    summary = {
        "isolates": analyze_isolates(),
        "patients": analyze_patients(),
    }
    with open(PROC / "spidaar_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)
    print("SPIDAAR isolates n =", summary["isolates"]["n"], "| patients n =", summary["patients"]["n"])
    print(f"Saved summaries and figures to {PROC} and {FIG}")


if __name__ == "__main__":
    main()
