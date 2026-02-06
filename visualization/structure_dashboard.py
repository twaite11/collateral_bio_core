"""
Generate Structure Dashboard HTML: 3D viewer + domain coloring + homology table.
Reads: input FASTA, InterPro domains.json, homology_scores.json, ColabFold PDB paths.
Serve the HTML via HTTP (e.g. python -m http.server) to load PDBs.
"""
import argparse
import json
import re
from pathlib import Path

from Bio import SeqIO

HEPN_REGEX = re.compile(r"R.{4,6}H")


def count_hepn(seq_str: str) -> int:
    return len(HEPN_REGEX.findall(str(seq_str)))


def parse_family_id(seq_id: str) -> str:
    parts = str(seq_id).split("_")
    return parts[0] if parts else seq_id


def find_pdb_paths(structures_dir: Path) -> dict:
    """Scan ColabFold output for {seq_id: relative_path}."""
    result = {}
    for sub in structures_dir.iterdir():
        if not sub.is_dir():
            continue
        for pdb in sub.glob("*_rank_001*.pdb"):
            name = pdb.stem.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
            rel = pdb.relative_to(structures_dir.parent.parent.parent)  # from project root
            result[name] = str(Path("data") / rel).replace("\\", "/")
    for pdb in structures_dir.glob("*.pdb"):
        result[pdb.stem] = "data/structure_pipeline/structures/colabfold/" + pdb.name
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate structure dashboard with 3D viewer and homology table."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="Input FASTA (2-3 HEPN)",
    )
    parser.add_argument(
        "--domains",
        default="data/structure_pipeline/interpro/domains.json",
        help="InterPro domains JSON",
    )
    parser.add_argument(
        "--homology",
        default="data/structure_pipeline/homology_scores.json",
        help="Homology scores JSON",
    )
    parser.add_argument(
        "--structures-dir",
        default="data/structure_pipeline/structures/colabfold",
        help="ColabFold output directory",
    )
    parser.add_argument(
        "--output",
        default="visualization/structure_dashboard.html",
        help="Output HTML path",
    )
    args = parser.parse_args()

    proj = Path(__file__).resolve().parent.parent
    input_path = proj / args.input
    domains_path = proj / args.domains
    homology_path = proj / args.homology
    struct_dir = proj / args.structures_dir
    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = proj / out_path

    # Load sequences
    seqs = {}
    if input_path.exists():
        for rec in SeqIO.parse(str(input_path), "fasta"):
            seqs[rec.id] = {
                "length": len(rec.seq),
                "hepn_count": count_hepn(str(rec.seq)),
                "family": parse_family_id(rec.id),
            }
    else:
        seqs = {}

    # Load domains
    domains = {}
    if domains_path.exists():
        with open(domains_path) as f:
            domains = json.load(f)

    # Load homology
    homology = {}
    if homology_path.exists():
        with open(homology_path) as f:
            homology = json.load(f)

    # PDB paths (relative for fetch from project root)
    pdb_paths = find_pdb_paths(struct_dir) if struct_dir.exists() else {}

    # Build rows for table
    all_ids = sorted(set(seqs.keys()) | set(domains.keys()) | set(homology.keys()) | set(pdb_paths.keys()))
    rows = []
    for seq_id in all_ids:
        s = seqs.get(seq_id, {})
        d = domains.get(seq_id, [])
        h = homology.get(seq_id, {})
        cas13a = h.get("cas13a")
        cas13b = h.get("cas13b")
        best = None
        best_name = None
        if cas13a is not None and (best is None or (cas13a or 0) > best):
            best = cas13a
            best_name = "Cas13a"
        if cas13b is not None and (best is None or (cas13b or 0) > best):
            best = cas13b
            best_name = "Cas13b"
        rows.append({
            "id": seq_id,
            "length": s.get("length", 0),
            "hepn_count": s.get("hepn_count", 0),
            "family": s.get("family", parse_family_id(seq_id)),
            "domains": d,
            "cas13a": round(cas13a * 100, 1) if cas13a is not None else None,
            "cas13b": round(cas13b * 100, 1) if cas13b is not None else None,
            "best_match": best_name,
            "best_pct": round((best or 0) * 100, 1) if best is not None else None,
            "pdb_path": pdb_paths.get(seq_id),
        })

    # Sort by best homology desc
    rows.sort(key=lambda r: (r["best_pct"] or 0, r["id"]), reverse=True)

    # Domain color map for 3Dmol
    domain_colors = {
        "PF05168": "#22d3ee",   # HEPN
        "PF01228": "#f97316",   # HEL
        "PF18456": "#94a3b8",   # WYL
    }

    data_json = json.dumps({
        "rows": rows,
        "domain_colors": domain_colors,
    })

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Structure Dashboard - 2-3 HEPN Cas13</title>
  <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    :root {{ --bg:#0f172a; --card:#1e3a5f; --accent:#38bdf8; --hepn:#22d3ee; --text:#f8fafc; --muted:#94a3b8; }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; padding: 1.5rem; font-family: 'DM Sans', system-ui; background: var(--bg); color: var(--text); min-height: 100vh; }}
    h1 {{ font-size: 1.5rem; margin: 0 0 0.5rem; }}
    .sub {{ color: var(--muted); font-size: 0.9rem; margin-bottom: 1rem; }}
    .layout {{ display: grid; grid-template-columns: 1fr 400px; gap: 1.5rem; }}
    @media (max-width: 900px) {{ .layout {{ grid-template-columns: 1fr; }} }}
    .viewer {{ width: 100%; height: 480px; background: #1e293b; border-radius: 8px; }}
    .table-wrap {{ overflow-x: auto; max-height: 70vh; overflow-y: auto; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
    th, td {{ padding: 0.4rem 0.6rem; text-align: left; border-bottom: 1px solid #334155; }}
    th {{ color: var(--muted); font-weight: 500; position: sticky; top: 0; background: var(--bg); }}
    tr:hover {{ background: #1e3a5f; }}
    tr.selected {{ background: #1e3a5f; outline: 1px solid var(--accent); }}
    .pct {{ font-variant-numeric: tabular-nums; }}
    .no-pdb {{ color: var(--muted); font-style: italic; }}
    .load-hint {{ color: var(--muted); font-size: 0.8rem; margin-top: 0.5rem; }}
  </style>
</head>
<body>
  <h1>Structure Dashboard – 2-3 HEPN Cas13</h1>
  <p class="sub">Folding homology vs Cas13a/13b · Domain-colored 3D structure · Serve via HTTP to load PDBs</p>

  <div class="layout">
    <div>
      <div id="molviewer" class="viewer"></div>
      <p class="load-hint">Select a row to load structure. Serve from project root: python -m http.server 8000</p>
      <div id="domain-legend" style="margin-top:0.5rem;font-size:0.85rem;color:var(--muted)"></div>
    </div>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>ID</th>
            <th>Family</th>
            <th>Len</th>
            <th>HEPN</th>
            <th>Cas13a %</th>
            <th>Cas13b %</th>
            <th>Best</th>
          </tr>
        </thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </div>

  <script>
    const DATA = {data_json};

    function renderTable() {{
      const tbody = document.getElementById("tbody");
      tbody.innerHTML = "";
      DATA.rows.forEach((r, i) => {{
        const tr = document.createElement("tr");
        tr.dataset.id = r.id;
        tr.dataset.index = i;
        tr.innerHTML = `<td>${{r.id}}</td><td>${{r.family}}</td><td>${{r.length}}</td><td>${{r.hepn_count}}</td>
          <td class="pct">${{r.cas13a != null ? r.cas13a + '%' : '-'}}</td>
          <td class="pct">${{r.cas13b != null ? r.cas13b + '%' : '-'}}</td>
          <td>${{r.best_match || '-'}}</td>`;
        tr.addEventListener("click", () => selectRow(i));
        tbody.appendChild(tr);
      }});
    }}

    function selectRow(index) {{
      document.querySelectorAll("#tbody tr").forEach(tr => tr.classList.remove("selected"));
      const tr = document.querySelector(`#tbody tr[data-index="${{index}}"]`);
      if (tr) tr.classList.add("selected");
      loadStructure(DATA.rows[index]);
    }}

    function loadStructure(row) {{
      const viewer = document.getElementById("molviewer");
      const div = document.createElement("div");
      div.className = "viewer";
      div.style.width = "100%";
      div.style.height = "480px";
      viewer.parentNode.replaceChild(div, viewer);
      div.id = "molviewer";

      const view = $3Dmol.createViewer(div, {{ backgroundColor: "#1e293b" }});

      if (!row.pdb_path) {{
        view.addLabel("No PDB available", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
        return;
      }}

      // Fetch PDB - path relative to project root when served via http.server
      const base = window.location.pathname.includes("/visualization/") ? "../" : "./";
      const url = base + (row.pdb_path || "");
      fetch(url).then(res => res.text()).then(pdb => {{
        view.addModel(pdb, "pdb");
        view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});

        // Domain coloring from InterPro
        (row.domains || []).forEach(d => {{
          const sig = (d.signature || "").split(".")[0];
          const color = DATA.domain_colors[sig] || "#64748b";
          view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: color }} }});
        }});

        view.zoomTo();
        view.render();

        // Legend
        const legend = document.getElementById("domain-legend");
        const legMap = {{ "PF05168": "HEPN", "PF01228": "HEL", "PF18456": "WYL" }};
        const used = [...new Set((row.domains || []).map(d => (d.signature || "").split(".")[0]))];
        legend.innerHTML = used.map(s => `<span style="color:${{DATA.domain_colors[s] || '#64748b'}}">&#9679;</span> ${{legMap[s] || s}} `).join(" ") || "";
      }}).catch(() => {{
        view.addLabel("Failed to load PDB", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
      }});
    }}

    renderTable();
    if (DATA.rows.length) selectRow(0);
  </script>
</body>
</html>
'''

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[SUCCESS] Wrote {out_path}")
    print(f"   [*] Serve via: python -m http.server 8000  (from project root)")
    print(f"   [*] Open: http://localhost:8000/visualization/structure_dashboard.html")


if __name__ == "__main__":
    main()
