"""Build focused post-npj revision figures from audited project tables.

The script does not refit the original dependency or proteogenomic analyses.
It redraws their frozen summaries and replots published source values from
Ye et al. (Nature Communications, 2026) without new significance testing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import openpyxl
import pandas as pd


INK = "#17324D"
BLUE = "#22577A"
TEAL = "#2A9D8F"
GOLD = "#E9C46A"
ORANGE = "#E76F51"
PURPLE = "#6D597A"
GREY = "#657786"
LIGHT = "#EAF0F4"
RED = "#B23A48"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig: plt.Figure, stem: Path) -> dict:
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=450, facecolor="white")
    fig.savefig(stem.parent / f"{stem.name}_preview.png", dpi=170, facecolor="white")
    plt.close(fig)
    doc = fitz.open(stem.with_suffix(".pdf"))
    page = doc[0]
    spans = [
        span
        for block in page.get_text("dict")["blocks"]
        if block["type"] == 0
        for line in block["lines"]
        for span in line["spans"]
    ]
    check = {
        "pdf": stem.with_suffix(".pdf").name,
        "width_mm": page.rect.width * 25.4 / 72,
        "height_mm": page.rect.height * 25.4 / 72,
        "minimum_font_pt": min(span["size"] for span in spans),
        "text_outside_page": [
            span["text"]
            for span in spans
            if not page.rect.contains(fitz.Rect(span["bbox"]))
        ],
        "embedded_rasters": len(page.get_images(full=True)),
    }
    doc.close()
    return check


def panel_label(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(-0.08, 1.07, label, transform=ax.transAxes, fontsize=11, weight="bold")
    ax.set_title(title, loc="left", fontsize=9.2, weight="bold", pad=10)


def rounded_box(ax, x, y, w, h, title, body, color, title_fs=8.2, body_fs=7.2):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.01,rounding_size=0.02",
            facecolor="white",
            edgecolor=color,
            linewidth=1.25,
        )
    )
    ax.text(x + w / 2, y + h * 0.73, title, ha="center", va="center", weight="bold", color=color, fontsize=title_fs)
    ax.text(x + w / 2, y + h * 0.36, body, ha="center", va="center", color=INK, fontsize=body_fs, linespacing=1.35)


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=10,
            color=GREY,
            linewidth=1.1,
        )
    )


def figure1(output_dir: Path) -> dict:
    fig, ax = plt.subplots(figsize=(180 / 25.4, 118 / 25.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.96, "Focused biological question", ha="center", va="center", fontsize=11, weight="bold", color=INK)
    ax.text(
        0.5,
        0.895,
        "Does a structure-selected cBAF panel show reproducible residual protein co-variation in lung cancer?",
        ha="center",
        va="center",
        fontsize=8.4,
        color=GREY,
    )
    xs = [0.035, 0.275, 0.515, 0.755]
    titles = ["Define", "Test", "Replicate", "Interpret"]
    bodies = [
        "6LTJ/6LTH contacts\nSix selected proteins\nBinary contact sensitivity",
        "Lung CRISPR co-dependency\n694 matched cell lines\nRNA/CNV/lineage adjustment",
        "CPTAC LUAD (n=108)\nCPTAC LUSC (n=107)\nProteome-wide null panels",
        "Published ARID1A knockout\nH1299/A549 IDR rescue\nARID1B fitness boundary",
    ]
    colors = [PURPLE, BLUE, TEAL, ORANGE]
    for x, title, body, color in zip(xs, titles, bodies, colors):
        rounded_box(ax, x, 0.50, 0.205, 0.27, title, body, color)
    for x1, x2 in zip(xs[:-1], xs[1:]):
        arrow(ax, x1 + 0.205, 0.635, x2, 0.635)
    rounded_box(
        ax,
        0.035,
        0.14,
        0.445,
        0.23,
        "ARID1A state is not interchangeable",
        "Alteration: any annotated variant\nPredicted loss of function: truncating/splice variant\nEngineered loss: verified knockout\nProtein-low: abundance proxy only",
        RED,
    )
    rounded_box(
        ax,
        0.52,
        0.14,
        0.445,
        0.23,
        "Inference boundary",
        "Residual co-variation is not a direct-contact effect\nCo-dependency is not loss-conditioned lethality\nARID1B support does not imply whole-cBAF dependence\nNo therapeutic marker is claimed",
        GREY,
    )
    ax.text(
        0.5,
        0.045,
        "Primary endpoint: selected-panel minus broad-BAF-background median correlation of linear-model residuals.",
        ha="center",
        fontsize=7.2,
        color=GREY,
    )
    fig.subplots_adjust(0.015, 0.02, 0.985, 0.98)
    return save(fig, output_dir / "Figure1_focused")


def cryoem_crop(pdf_path: Path) -> np.ndarray:
    doc = fitz.open(pdf_path)
    page = doc[0]
    clip = fitz.Rect(0, 0, page.rect.width, page.rect.height * 0.30)
    pix = page.get_pixmap(matrix=fitz.Matrix(3.2, 3.2), clip=clip, alpha=False)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    doc.close()
    return arr[:, :, :3]


def figure2(project_root: Path, revision_dir: Path, output_dir: Path) -> dict:
    contacts_file = project_root / "results_tables" / "structure_validation_r5" / "structure_contact_cutoff_sensitivity.csv"
    global_file = project_root / "results_tables" / "structure_validation_r5" / "structure_global_validation.csv"
    groups_file = project_root / "results_tables" / "structure_validation_r5" / "structural_group_sensitivity.csv"
    contacts = pd.read_csv(contacts_file)
    global_q = pd.read_csv(global_file).set_index("pdbid")
    groups = pd.read_csv(groups_file)
    genes = ["SMARCC2", "SMARCA4", "SMARCD1", "SMARCB1", "SMARCE1", "DPF2"]
    observed = contacts[contacts["gene"].isin(genes)].copy()
    assert len(observed) == 36
    assert observed["contact_call_after_Q_filter"].eq(1).all()

    fig = plt.figure(figsize=(180 / 25.4, 155 / 25.4))
    ax_a = fig.add_axes([0.07, 0.71, 0.40, 0.20])
    panel_label(ax_a, "a", "ARID1A molecular regions")
    ax_a.set_xlim(0, 2285)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    regions = [
        (1, 1012, "IDR1-3", LIGHT),
        (1013, 1110, "ARID", GOLD),
        (1111, 1634, "cIDR", TEAL),
        (1635, 2285, "structured GR", PURPLE),
    ]
    for start, end, label, color in regions:
        ax_a.add_patch(Rectangle((start, 0.42), end - start + 1, 0.24, facecolor=color, edgecolor="white"))
        ax_a.text((start + end) / 2, 0.54, label, ha="center", va="center", fontsize=7.2, weight="bold")
        boundary_y = 0.29 if start in {1, 1013, 1635} else 0.18
        ax_a.text(start, boundary_y, str(start), ha="center", fontsize=6.4, color=GREY)
    ax_a.text(2285, 0.30, "2285", ha="center", fontsize=6.7, color=GREY)
    ax_a.text(0.00, 0.02, "IDRs: condensation / TF recruitment", transform=ax_a.transAxes, fontsize=6.2, color=TEAL, ha="left")
    ax_a.text(0.98, 0.02, "GR: cBAF contacts", transform=ax_a.transAxes, fontsize=6.2, color=PURPLE, ha="right")

    ax_b = fig.add_axes([0.52, 0.64, 0.43, 0.31])
    panel_label(ax_b, "b", "Deposited density and fitted model")
    ax_b.imshow(cryoem_crop(revision_dir / "figures" / "Figure3.pdf"))
    ax_b.axis("off")

    ax_c = fig.add_axes([0.09, 0.15, 0.40, 0.38])
    panel_label(ax_c, "c", "Binary contacts are robust to deposit and cutoff")
    matrix = []
    columns = []
    for pdb in ["6LTJ", "6LTH"]:
        for cutoff in [4.0, 4.5, 5.0]:
            columns.append(f"{pdb}\n{cutoff:.1f} A")
            sub = observed[(observed["pdbid"] == pdb) & np.isclose(observed["cutoff_A"], cutoff)].set_index("gene")
            matrix.append(sub.loc[genes, "contact_call_after_Q_filter"].to_numpy())
    mat = np.vstack(matrix).T
    ax_c.imshow(mat, cmap=matplotlib.colors.ListedColormap(["white", TEAL]), vmin=0, vmax=1, aspect="auto")
    ax_c.set_xticks(range(len(columns)), columns, fontsize=6.8)
    ax_c.set_yticks(range(len(genes)), genes, fontsize=7.2)
    for i in range(len(genes)):
        for j in range(len(columns)):
            ax_c.text(j, i, "yes", ha="center", va="center", fontsize=6.1, color="white", weight="bold")
    ax_c.set_xlabel("Contact retained after bottom-decile residue-Q removal", fontsize=6.8)
    ax_c.tick_params(length=0)
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    ax_d = fig.add_axes([0.62, 0.17, 0.32, 0.34])
    panel_label(ax_d, "d", "Structure-selected panel contrast")
    selected = groups[groups["analysis"] == "experimental_6LTJ_6LTH_contacts_only_6v8"].set_index("cohort")
    cohorts = ["lung", "pancancer"]
    labels = ["Lung\n(128 models)", "Pan-cancer\n(930 models)"]
    for idx, (cohort, label, color) in enumerate(zip(cohorts, labels, [BLUE, ORANGE])):
        row = selected.loc[cohort]
        value = row["median_difference"]
        low = row["bootstrap_95ci_low"]
        high = row["bootstrap_95ci_high"]
        ax_d.errorbar(value, idx, xerr=[[value - low], [high - value]], fmt="o", color=color, capsize=3, markersize=6)
    ax_d.axvline(0, color=GREY, linestyle=":", linewidth=1)
    ax_d.set_yticks([0, 1], labels)
    ax_d.invert_yaxis()
    ax_d.set_xlabel("Median Spearman rho difference\n(6 selected proteins minus 8 comparisons)", fontsize=7.2)
    ax_d.text(
        0.98,
        0.50,
        "Exact two-sided P: 0.0080 / 0.0127",
        transform=ax_d.transAxes,
        fontsize=6.6,
        color=GREY,
        ha="right",
        va="center",
    )
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.tick_params(labelsize=7)

    fig.text(
        0.5,
        0.035,
        f"All six calls persisted in 6LTJ/6LTH at 4.0-5.0 A; global Q-scores were {global_q.loc['6LTJ','global_Q_score']:.3f} and {global_q.loc['6LTH','global_Q_score']:.3f}.\n"
        "This exploratory contrast does not isolate direct contact; SMARCC1 is not in the selected panel.",
        ha="center",
        fontsize=7.1,
        color=GREY,
    )
    return save(fig, output_dir / "Figure2_focused")


def read_block(sheet, start_row: int, end_row: int, columns: list[int]) -> np.ndarray:
    values = []
    for row in range(start_row, end_row + 1):
        values.append([float(sheet.cell(row=row, column=column).value) for column in columns])
    return np.asarray(values)


def molecular_source_values(source_xlsx: Path) -> tuple[pd.DataFrame, dict]:
    workbook = openpyxl.load_workbook(source_xlsx, read_only=True, data_only=True)
    h1299 = workbook["Figure 7"]
    a549 = workbook["Supplementary Figure 13"]
    # Published 96-hour values: H1299 Figure 7f rows 58-62, columns N-S;
    # A549 Supplementary Figure 13b rows 26-30, columns F-K.
    h_values = read_block(h1299, 58, 62, [14, 15, 16, 17, 18, 19])
    a_values = read_block(a549, 26, 30, [6, 7, 8, 9, 10, 11])
    groups = ["ARID1A", "Delta-cIDR", "Y-to-A", "GATA2+ARID1A", "GATA2+Delta-cIDR", "GATA2+Y-to-A"]
    rows = []
    for cell_line, matrix in [("H1299", h_values), ("A549", a_values)]:
        for col, group in enumerate(groups):
            for replicate, value in enumerate(matrix[:, col], start=1):
                rows.append({"cell_line": cell_line, "condition": group, "replicate": replicate, "relative_proliferation_96h": value})
    frame = pd.DataFrame(rows)
    summary = (
        frame.groupby(["cell_line", "condition"], sort=False)["relative_proliferation_96h"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    descriptive = {}
    for cell_line in ["H1299", "A549"]:
        means = summary[summary["cell_line"] == cell_line].set_index("condition")["mean"]
        descriptive[cell_line] = {
            "GATA2_gain_with_ARID1A": float(means["GATA2+ARID1A"] - means["ARID1A"]),
            "GATA2_gain_with_Delta_cIDR": float(means["GATA2+Delta-cIDR"] - means["Delta-cIDR"]),
            "GATA2_gain_with_Y_to_A": float(means["GATA2+Y-to-A"] - means["Y-to-A"]),
        }
    return frame, {"summary": summary.to_dict(orient="records"), "descriptive_GATA2_gains": descriptive}


def figure5(source_xlsx: Path, output_dir: Path) -> tuple[dict, dict]:
    frame, analysis = molecular_source_values(source_xlsx)
    conditions = ["ARID1A", "Delta-cIDR", "Y-to-A", "GATA2+ARID1A", "GATA2+Delta-cIDR", "GATA2+Y-to-A"]
    labels = ["WT", "Delta-cIDR", "Y-to-A", "GATA2+WT", "GATA2+Delta-cIDR", "GATA2+Y-to-A"]
    colors = [BLUE, TEAL, PURPLE, BLUE, TEAL, PURPLE]
    markers = ["o", "o", "o", "D", "D", "D"]
    fig = plt.figure(figsize=(180 / 25.4, 125 / 25.4))
    for panel, cell_line, left in [("a", "H1299", 0.08), ("b", "A549", 0.56)]:
        ax = fig.add_axes([left, 0.43, 0.37, 0.44])
        ax.text(-0.08, 1.06, panel, transform=ax.transAxes, fontsize=11, weight="bold")
        ax.set_title(f"{cell_line} rescue (96 h)", loc="left", fontsize=9.0, weight="bold", pad=9)
        subset = frame[frame["cell_line"] == cell_line]
        for idx, (condition, color, marker) in enumerate(zip(conditions, colors, markers)):
            values = subset[subset["condition"] == condition]["relative_proliferation_96h"].to_numpy()
            jitter = np.linspace(-0.08, 0.08, len(values))
            ax.scatter(idx + jitter, values, s=22, marker=marker, facecolor=color if marker == "D" else "white", edgecolor=color, linewidth=1)
            ax.errorbar(idx, values.mean(), yerr=values.std(ddof=1), fmt="_", color=INK, markersize=13, capsize=2, linewidth=1)
        ax.set_xticks(range(len(labels)), labels, rotation=35, ha="right")
        ax.set_ylabel("Relative proliferation")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=6.8)
        ax.axvline(2.5, color=LIGHT, linewidth=1.2)
        ax.text(0.20, 0.96, "without GATA2", transform=ax.transAxes, fontsize=6.8, color=GREY, ha="center")
        ax.text(0.76, 0.96, "with GATA2", transform=ax.transAxes, fontsize=6.8, color=GREY, ha="center")
    ax_c = fig.add_axes([0.08, 0.075, 0.85, 0.22])
    ax_c.axis("off")
    boxes = [
        (0.00, "Molecular anchor", "cIDR-dependent TF recruitment\nin published wet-lab studies", TEAL),
        (0.35, "Primary observation", "Six structure-selected proteins show\nresidual co-variation in LUAD and LUSC", BLUE),
        (0.70, "Functional boundary", "Published knockout/rescue gives context;\nwhole-module dependence is unproven", ORANGE),
    ]
    for x, title, body, color in boxes:
        rounded_box(ax_c, x, 0.05, 0.28, 0.82, title, body, color, title_fs=7.6, body_fs=6.5)
    arrow(ax_c, 0.28, 0.46, 0.35, 0.46)
    arrow(ax_c, 0.63, 0.46, 0.70, 0.46)
    fig.text(
        0.5,
        0.015,
        "All five published replicate values are shown; bars are mean +/- SD. No new significance test was performed.",
        ha="center",
        fontsize=7.0,
        color=GREY,
    )
    check = save(fig, output_dir / "Figure5_molecular_validation")
    frame.to_csv(output_dir / "Figure5_molecular_validation_source.csv", index=False)
    (output_dir / "Figure5_molecular_validation_summary.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    return check, analysis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--source-xlsx", type=Path, required=True)
    parser.add_argument("--revision-dir", type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.revision_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    checks = [figure1(output_dir), figure2(args.project_root, args.revision_dir, output_dir)]
    check5, analysis = figure5(args.source_xlsx, output_dir)
    checks.append(check5)
    manifest = {
        "script_sha256": sha256(Path(__file__)),
        "source_xlsx": str(args.source_xlsx),
        "source_xlsx_sha256": sha256(args.source_xlsx),
        "frozen_project_sources": {
            "structure_contact_cutoff_sensitivity.csv": sha256(args.project_root / "results_tables" / "structure_validation_r5" / "structure_contact_cutoff_sensitivity.csv"),
            "structure_global_validation.csv": sha256(args.project_root / "results_tables" / "structure_validation_r5" / "structure_global_validation.csv"),
            "structural_group_sensitivity.csv": sha256(args.project_root / "results_tables" / "structure_validation_r5" / "structural_group_sensitivity.csv"),
        },
        "new_significance_tests": False,
        "checks": checks,
        "external_molecular_descriptive_summary": analysis,
    }
    (args.revision_dir / "focused_figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
