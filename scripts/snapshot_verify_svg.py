#!/usr/bin/env python3
"""Render a verify.json into a self-contained SVG snapshot of the verification matrix.

Used to produce the README preview image (docs/img/verify_matrix.svg). SVG (not a PNG
screenshot) so it needs no headless browser, stays crisp, and renders on GitHub. Colours
mirror the HTML report: green on-target, red off-target, grey dropout/uncertain.

Usage:  python scripts/snapshot_verify_svg.py <verify.json> <out.svg>
"""
from __future__ import annotations

import json
import sys
from html import escape

GREEN, RED, GREY, FG, LINE, HEAD_BG = "#1a7f37", "#cf222e", "#6e7781", "#1f2328", "#d0d7de", "#f6f8fa"
CW, PAD, ROW_H, HEAD_H = 7.8, 12.0, 30.0, 34.0  # monospace char width, cell pad, row/header height
FONT = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"


def cell_segments(c: dict) -> list[tuple[str, str, bool]]:
    """(text, colour, italic) segments for one matrix cell."""
    if c["status"] == "uncertain":
        return [("?", GREY, False)]
    # Before the empty-product test: a cell capped by search.max_binding_sites enumerates no
    # products, so it would otherwise render as "dropout" -- reading a promiscuous primer as
    # a failure to amplify, which is the exact mislabel the cap exists to prevent.
    if c.get("site_cap"):
        return [(f">{c['site_cap']} binding sites", RED, False)]
    if not c["on_target"] and not c["off_target"]:
        return [("dropout", GREY, True)]
    segs: list[tuple[str, str, bool]] = []
    if c["on_target"]:
        segs.append((",".join(map(str, c["on_target"])), GREEN, False))
    if c["off_target"]:
        if segs:
            segs.append(("  ", FG, False))
        segs.append(("off:" + ",".join(map(str, c["off_target"])), RED, False))
    return segs


def plain(segs) -> str:
    return "".join(t for t, _, _ in segs)


def build(data: dict) -> str:
    haps = data["haplotypes"]
    header = ["primer_id", "expected bp", *haps]
    # each body row: list of cells, each a list of (text, colour, italic) segments
    body: list[list[list[tuple[str, str, bool]]]] = []
    for r in data["rows"]:
        cm = {c["haplotype_id"]: c for c in r["cells"]}
        row = [[(r["primer_id"], FG, False)], [(str(r["expected_size"]), FG, False)]]
        row += [cell_segments(cm[h]) for h in haps]
        body.append(row)

    ncol = len(header)
    # column width = widest plain text in that column (header or any body cell)
    widths = []
    for j in range(ncol):
        longest = len(header[j])
        for row in body:
            longest = max(longest, len(plain(row[j])))
        widths.append(longest * CW + 2 * PAD)
    xs = [0.0]
    for w in widths:
        xs.append(xs[-1] + w)
    table_w = xs[-1]

    m = 20.0                       # outer margin
    title_y, legend_y = 34.0, 60.0
    top = 78.0                     # table top
    table_h = HEAD_H + ROW_H * len(body)
    W = table_w + 2 * m
    H = top + table_h + m

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
        f'viewBox="0 0 {W:.0f} {H:.0f}" font-family="{escape(FONT)}">',
        f'<rect width="{W:.0f}" height="{H:.0f}" fill="#ffffff"/>',
        f'<text x="{m}" y="{title_y}" font-size="19" font-weight="700" fill="{FG}">'
        f'Primer verification across the pangenome</text>',
    ]
    # legend
    lx = m
    for text, colour in [("correct size", GREEN), ("off-target", RED),
                         ("dropout", GREY), ("? not projectable", GREY)]:
        out.append(f'<text x="{lx:.1f}" y="{legend_y}" font-size="12.5" fill="{colour}" '
                   f'font-weight="600">{escape(text)}</text>')
        lx += (len(text) + 3) * 7.0

    gx = m                                    # grid origin x
    gy = top
    # header background
    out.append(f'<rect x="{gx:.1f}" y="{gy:.1f}" width="{table_w:.1f}" height="{HEAD_H:.1f}" '
               f'fill="{HEAD_BG}"/>')

    def text_cell(cx, cy, segs, weight="400"):
        # left-aligned monospace text with per-segment colour via tspans
        parts = [f'<text x="{cx + PAD:.1f}" y="{cy:.1f}" font-size="13" font-weight="{weight}">']
        for t, colour, italic in segs:
            style = ' font-style="italic"' if italic else ""
            # preserve spaces between segments
            parts.append(f'<tspan xml:space="preserve" fill="{colour}"{style}>{escape(t)}</tspan>')
        parts.append("</text>")
        return "".join(parts)

    # header row
    hy = gy + HEAD_H / 2 + 5
    for j, htext in enumerate(header):
        out.append(text_cell(gx + xs[j], hy, [(htext, FG, False)], weight="700"))
    # body rows
    for i, row in enumerate(body):
        ry = gy + HEAD_H + i * ROW_H
        cy = ry + ROW_H / 2 + 5
        for j, segs in enumerate(row):
            weight = "700" if j == 0 else "400"
            out.append(text_cell(gx + xs[j], cy, segs, weight=weight))

    # grid lines
    for j in range(ncol + 1):
        x = gx + xs[j]
        out.append(f'<line x1="{x:.1f}" y1="{gy:.1f}" x2="{x:.1f}" y2="{gy + table_h:.1f}" '
                   f'stroke="{LINE}"/>')
    for i in range(len(body) + 2):
        y = gy + (HEAD_H if i else 0) + max(0, i - 1) * ROW_H
        out.append(f'<line x1="{gx:.1f}" y1="{y:.1f}" x2="{gx + table_w:.1f}" y2="{y:.1f}" '
                   f'stroke="{LINE}"/>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    data = json.loads(open(sys.argv[1]).read())
    open(sys.argv[2], "w").write(build(data))
    print(f"wrote {sys.argv[2]}")


if __name__ == "__main__":
    main()
