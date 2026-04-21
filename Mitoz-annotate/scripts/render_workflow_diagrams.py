#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render bilingual workflow diagrams for batch_mitoz.py.

The SVG files are written as UTF-8, which prevents encoding mojibake.  They
also include a CJK font stack so Chinese labels can use common system Chinese
fonts without relying on matplotlib, reducing missing-glyph boxes.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "flowcharts"

# Chinese font support: SVG viewers will choose the first installed family.
# UTF-8 handles encoding; this font stack handles CJK glyph coverage when any
# common Chinese font is installed on Linux/Windows/macOS.
CJK_FONT_STACK = (
    '"Noto Sans CJK SC", "Source Han Sans SC", "Microsoft YaHei", '
    '"SimHei", "WenQuanYi Micro Hei", "PingFang SC", Arial, sans-serif'
)


Node = Tuple[str, Sequence[str]]


def _text_lines(lines: Sequence[str], x: int, y: int, size: int = 16) -> str:
    tspans: List[str] = []
    for i, line in enumerate(lines):
        dy = "0" if i == 0 else "1.35em"
        tspans.append(
            f'<tspan x="{x}" dy="{dy}">{html.escape(line)}</tspan>'
        )
    return (
        f'<text x="{x}" y="{y}" font-size="{size}" '
        f'text-anchor="middle" dominant-baseline="middle">'
        + "".join(tspans)
        + "</text>"
    )


def render_svg(title: str, subtitle: str, nodes: Iterable[Node], out_path: Path) -> None:
    nodes = list(nodes)
    width = 980
    node_x = 140
    node_w = 700
    node_h = 72
    gap = 28
    top = 112
    height = top + len(nodes) * node_h + (len(nodes) - 1) * gap + 72
    center_x = node_x + node_w // 2

    parts: List[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">',
        "<defs>",
        "<style>",
        f"text {{ font-family: {CJK_FONT_STACK}; fill: #17202a; }}",
        ".title { font-weight: 700; }",
        ".subtitle { fill: #566573; }",
        ".node { fill: #f8fbff; stroke: #2f6f9f; stroke-width: 2; rx: 8; }",
        ".node-accent { fill: #eef7f1; stroke: #2f7d50; stroke-width: 2; rx: 8; }",
        ".arrow { stroke: #586069; stroke-width: 2.4; fill: none; marker-end: url(#arrowhead); }",
        "</style>",
        '<marker id="arrowhead" markerWidth="10" markerHeight="10" refX="8" refY="3" '
        'orient="auto" markerUnits="strokeWidth">',
        '<path d="M0,0 L0,6 L9,3 z" fill="#586069" />',
        "</marker>",
        "</defs>",
        '<rect x="0" y="0" width="980" height="100%" fill="#ffffff" />',
        f'<text class="title" x="{width // 2}" y="42" font-size="26" '
        f'text-anchor="middle">{html.escape(title)}</text>',
        f'<text class="subtitle" x="{width // 2}" y="72" font-size="15" '
        f'text-anchor="middle">{html.escape(subtitle)}</text>',
    ]

    for i, (_node_id, lines) in enumerate(nodes):
        y = top + i * (node_h + gap)
        css = "node-accent" if i in (0, len(nodes) - 1) else "node"
        parts.append(
            f'<rect class="{css}" x="{node_x}" y="{y}" width="{node_w}" height="{node_h}" />'
        )
        parts.append(_text_lines(lines, center_x, y + node_h // 2 - 8))
        if i < len(nodes) - 1:
            y1 = y + node_h
            y2 = y + node_h + gap - 6
            parts.append(f'<path class="arrow" d="M {center_x} {y1 + 5} L {center_x} {y2}" />')

    parts.append("</svg>")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> int:
    zh_nodes: List[Node] = [
        ("input", ["输入数据", "FASTA / GenBank；FASTQ 双端可选"]),
        ("preflight", ["输入识别与预检查", "UTF-8、Biopython、MitoZ、Snakemake"]),
        ("gb_convert", ["GenBank 去重与转换", "保留 B64_2 与 B64 为不同编号"]),
        ("internal", ["生成内部 FASTA", "S00001...；记录 id_map.tsv"]),
        ("phase1", ["第一轮 MitoZ 注释", "输出 01.mitoz_anno1/ 和 sample.log"]),
        ("reorient", ["按目标基因重定向", "默认 cox1；保留 linear / circular"]),
        ("phase2", ["可选第二轮重注释", "if_reannotation_2=yes 时运行"]),
        ("finalize", ["最终 GenBank 生成", "ID 重写、metadata 恢复、/organism= 同步"]),
        ("summary", ["统计与问题清单", "annotation-summary.tsv、failed/mismatch 列表"]),
        ("merge", ["合并最终结果", "04.final_gb/*.gbf 与 all.final.gbf"]),
    ]

    en_nodes: List[Node] = [
        ("input", ["Input data", "FASTA / GenBank; optional paired FASTQ"]),
        ("preflight", ["Detect input and run preflight", "UTF-8, Biopython, MitoZ, Snakemake"]),
        ("gb_convert", ["Deduplicate and convert GenBank", "Keep local IDs such as B64_2 distinct from B64"]),
        ("internal", ["Create internal FASTA files", "S00001... plus id_map.tsv"]),
        ("phase1", ["MitoZ annotation, phase 1", "Writes 01.mitoz_anno1/ and sample.log files"]),
        ("reorient", ["Reorient by target gene", "Default cox1; preserve linear / circular topology"]),
        ("phase2", ["Optional phase-2 re-annotation", "Runs when if_reannotation_2=yes"]),
        ("finalize", ["Build final GenBank files", "Rewrite IDs, restore metadata, synchronize /organism="]),
        ("summary", ["Write statistics and problem lists", "annotation-summary.tsv and failed/mismatch files"]),
        ("merge", ["Merge final outputs", "04.final_gb/*.gbf and all.final.gbf"]),
    ]

    render_svg(
        "batch_mitoz.py 工作流程",
        "UTF-8 避免编码乱码；CJK 字体栈降低中文缺字方块风险",
        zh_nodes,
        OUT_DIR / "batch_mitoz_workflow.zh.svg",
    )
    render_svg(
        "batch_mitoz.py Workflow",
        "UTF-8 SVG with CJK-capable font stack for bilingual labels",
        en_nodes,
        OUT_DIR / "batch_mitoz_workflow.en.svg",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
