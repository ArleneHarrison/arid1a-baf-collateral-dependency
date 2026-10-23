# -*- coding: utf-8 -*-
"""Render the deposited EMD-0974 density beside its fitted 6LTJ model.

The density panel is a maximum-intensity projection of the deposited 3D map,
not a raw particle micrograph. Voxels below the depositor's primary contour
level (0.48) are masked; no smoothing, local contrast adjustment, denoising or
selective editing is applied. The companion coordinate render is retained from
the existing manuscript figure so that the experimental density and its fitted
model are visually distinct.
"""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import mrcfile
import numpy as np
from scipy import ndimage


ROOT = Path(__file__).resolve().parents[2]
MAP = ROOT / "data" / "raw" / "emdb" / "EMD-0974" / "emd_0974.map"
MAP_GZ = MAP.with_suffix(".map.gz")
MODEL_RENDER = ROOT / "results_figures" / "fig4_cryoem_white.png"
OUT = ROOT / "results_figures"
CONTOUR = 0.48
WIDTH_MM = 162


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def crop_bounds(mask: np.ndarray, margin_fraction: float = 0.06) -> tuple[int, int, int, int]:
    yy, xx = np.where(mask)
    if not len(xx):
        raise ValueError(f"No projected voxel reaches contour {CONTOUR}")
    height, width = mask.shape
    margin = max(8, int(max(xx.max() - xx.min(), yy.max() - yy.min()) * margin_fraction))
    x0 = max(0, int(xx.min()) - margin)
    x1 = min(width, int(xx.max()) + margin + 1)
    y0 = max(0, int(yy.min()) - margin)
    y1 = min(height, int(yy.max()) + margin + 1)
    return x0, x1, y0, y1


def main() -> None:
    if not MAP.is_file():
        raise FileNotFoundError(f"Decompress {MAP_GZ} to {MAP} before rendering")
    if not MODEL_RENDER.is_file():
        raise FileNotFoundError(MODEL_RENDER)

    with mrcfile.mmap(MAP, permissive=True, mode="r") as volume:
        # MRC arrays are indexed z, y, x; projection along z preserves x as the
        # horizontal axis, so the scale bar uses voxel_size.x.
        projection = np.nanmax(volume.data, axis=0).astype(np.float32)
        voxel_x = float(volume.voxel_size.x)
        shape_zyx = tuple(int(v) for v in volume.data.shape)

    contour_mask = projection >= CONTOUR
    labels, component_count = ndimage.label(contour_mask)
    component_sizes = np.bincount(labels.ravel())
    component_sizes[0] = 0
    largest_label = int(component_sizes.argmax())
    mask = labels == largest_label
    x0, x1, y0, y1 = crop_bounds(mask)
    projected = projection[y0:y1, x0:x1]
    retained = mask[y0:y1, x0:x1]
    shown = np.ma.masked_where(~retained, projected)
    vmax = float(np.nanpercentile(projected[retained], 99.5))
    if not np.isfinite(vmax) or vmax <= CONTOUR:
        raise ValueError("Invalid density display range")

    model = mpimg.imread(MODEL_RENDER)
    mh, mw = model.shape[:2]
    # The legacy render contains a baked-in title and footer.  Reuse its
    # unchanged scientific content, but separate the legend from the molecular
    # render so that both main panels can occupy equal-width presentation boxes.
    # This is a layout-only crop; no density or structure pixels are adjusted.
    model_legend = model[
        int(mh * 0.085):int(mh * 0.345),
        int(mw * 0.015):int(mw * 0.305),
    ]
    model_structure = model[
        int(mh * 0.185):int(mh * 0.895),
        int(mw * 0.205):int(mw * 0.805),
    ]

    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    width_in = WIDTH_MM / 25.4
    fig = plt.figure(figsize=(width_in, width_in * 0.46), facecolor="white")

    # Manual equal-width panel boxes prevent off-axis notes from moving either
    # title.  Panel letters and titles are figure-level text at one shared y.
    title_y = 0.955
    fig.text(0.025, title_y, "a", ha="left", va="top", fontsize=8.6,
             fontweight="bold")
    fig.text(0.060, title_y, "Deposited cryo-EM density", ha="left", va="top",
             fontsize=8.2, fontweight="bold")
    fig.text(0.515, title_y, "b", ha="left", va="top", fontsize=8.6,
             fontweight="bold")
    fig.text(0.550, title_y, "Fitted coordinate model (PDB 6LTJ)", ha="left",
             va="top", fontsize=8.2, fontweight="bold")

    ax_density = fig.add_axes([0.035, 0.175, 0.430, 0.690])
    cmap = mpl.colormaps["Greys"].copy()
    cmap.set_bad("white")
    ax_density.imshow(shown, cmap=cmap, vmin=CONTOUR, vmax=vmax,
                      origin="lower", interpolation="nearest")
    ax_density.contour(retained.astype(float), levels=[0.5], colors="#4A4A4A",
                       linewidths=0.55)
    ax_density.axis("off")

    bar_angstrom = 100.0
    bar_pixels = bar_angstrom / voxel_x
    y_bar = shown.shape[0] * 0.075
    x_bar = shown.shape[1] * 0.08
    ax_density.plot([x_bar, x_bar + bar_pixels], [y_bar, y_bar],
                    color="black", linewidth=2.0, solid_capstyle="butt")
    ax_density.text(x_bar + bar_pixels / 2, y_bar + shown.shape[0] * 0.035,
                    "100 Å", ha="center", va="bottom", fontsize=6.5, color="black")
    fig.text(0.250, 0.075,
             "EMD-0974 · 3.7 Å · author contour 0.48\n"
             "largest connected density; no smoothing",
             ha="center", va="center", fontsize=6.2, color="#3F3F3F")

    ax_model = fig.add_axes([0.515, 0.135, 0.465, 0.730])
    ax_model.imshow(model_structure, interpolation="lanczos")
    ax_model.axis("off")

    # Compact legend occupies the naturally empty upper-left region of the
    # coordinate-model panel and does not determine the panel's data extent.
    ax_legend = fig.add_axes([0.525, 0.555, 0.155, 0.275], facecolor="white")
    ax_legend.imshow(model_legend, interpolation="lanczos")
    ax_legend.axis("off")

    png = OUT / "fig4_cryoem_density_composite.png"
    pdf = OUT / "fig4_cryoem_density_composite.pdf"
    fig.savefig(png, dpi=600, bbox_inches=None, facecolor="white")
    fig.savefig(pdf, dpi=600, bbox_inches=None, facecolor="white",
                metadata={"Title": "EMD-0974 cryo-EM density and fitted 6LTJ model"})
    plt.close(fig)

    manifest = {
        "figure_component": "Figure 4I-a",
        "source": {
            "accession": "EMD-0974",
            "map_url": (
                "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
                "EMD-0974/map/emd_0974.map.gz"
            ),
            "metadata_url": (
                "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
                "EMD-0974/header/emd-0974.xml"
            ),
            "accessed": date.today().isoformat(),
            "map_gz_sha256": sha256(MAP_GZ),
            "map_shape_zyx": shape_zyx,
            "voxel_size_x_angstrom": voxel_x,
        },
        "transformation": {
            "projection": "maximum intensity along MRC z axis",
            "mask_below_depositor_primary_contour": CONTOUR,
            "disconnected_density_policy": (
                "retain the largest 2D connected component of the projected "
                "author-contour mask; exclude disconnected dust islands"
            ),
            "projected_component_count": int(component_count),
            "retained_component_area_pixels": int(component_sizes[largest_label]),
            "upper_display_limit": "99.5th percentile of unmasked projected density",
            "interpolation": "nearest",
            "smoothing": "none",
            "local_adjustment": "none",
            "crop": {"x0": x0, "x1": x1, "y0": y0, "y1": y1,
                     "margin_fraction": 0.06},
            "layout": {
                "panel_width_fraction": 0.43,
                "shared_title_y_fraction": title_y,
                "model_render_policy": (
                    "unchanged legend and molecular-render pixels separated "
                    "from the legacy raster and repositioned; no resampling "
                    "interpreted as added structural detail"
                ),
                "model_legend_crop_fraction_xy": [0.015, 0.305, 0.085, 0.345],
                "model_structure_crop_fraction_xy": [0.205, 0.805, 0.185, 0.895],
            },
        },
        "outputs": [str(png.relative_to(ROOT)), str(pdf.relative_to(ROOT))],
    }
    (OUT / "fig4_cryoem_density_composite_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
