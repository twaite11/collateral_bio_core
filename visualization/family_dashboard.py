"""
Family Dashboard: Generates an HTML dashboard from fam_fasta.fasta
showing families, HEPN domain positions, and protein structure.
"""
import argparse
import json
import re
from pathlib import Path

from Bio import SeqIO

HEPN_REGEX = re.compile(r"R.{4,6}H")


def parse_family_id(seq_id: str) -> str:
    """Extract family ID from sequence ID (e.g. SN01_001 -> SN01)."""
    parts = str(seq_id).split("_")
    if len(parts) >= 1:
        return parts[0]
    return seq_id


def find_hepn_regions(seq_str: str) -> list[dict]:
    """Find all HEPN motifs; return list of {start, end, length, motif}."""
    regions = []
    for m in HEPN_REGEX.finditer(seq_str):
        start, end = m.start(), m.end()
        regions.append({
            "start": start,
            "end": end,
            "length": end - start,
            "motif": m.group(),
        })
    return regions


def load_families(fasta_path: str) -> dict:
    """Load FASTA, group by family, detect HEPN. Returns {family_id: [(seq_id, seq_len, hepn_regions), ...]}."""
    families = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(rec.seq)
        fam_id = parse_family_id(rec.id)
        hepn_regions = find_hepn_regions(seq_str)
        if fam_id not in families:
            families[fam_id] = []
        families[fam_id].append((rec.id, len(seq_str), hepn_regions))
    for fam_id in families:
        families[fam_id].sort(key=lambda x: x[0])
    return families


def generate_html(families: dict, output_path: str) -> None:
    """Generate self-contained HTML dashboard."""
    family_ids = sorted(families.keys(), key=lambda x: (len(x), x))
    total_members = sum(len(v) for v in families.values())

    # Build JSON for JS
    data = {}
    for fid in family_ids:
        members = []
        for seq_id, seq_len, hepn_regions in families[fid]:
            members.append({
                "id": seq_id,
                "length": seq_len,
                "hepn": hepn_regions,
            })
        data[fid] = members

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Family Dashboard - Deep Hits</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
  <style>
    :root {{
      --bg: #0f172a;
      --card: #1e3a5f;
      --card-hover: #334155;
      --accent: #38bdf8;
      --accent-dim: #0ea5e9;
      --hepn: #22d3ee;
      --text: #f8fafc;
      --text-muted: #94a3b8;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      padding: 2rem;
      font-family: 'DM Sans', system-ui, sans-serif;
      background: var(--bg);
      color: var(--text);
      min-height: 100vh;
    }}
    .header {{
      margin-bottom: 2rem;
    }}
    h1 {{
      font-size: 1.75rem;
      font-weight: 700;
      margin: 0 0 0.25rem 0;
      color: var(--text);
    }}
    .subtitle {{
      color: var(--text-muted);
      font-size: 0.95rem;
    }}
    .controls {{
      margin-bottom: 2rem;
    }}
    select {{
      padding: 0.6rem 1rem;
      font-size: 1rem;
      font-family: inherit;
      background: var(--card);
      color: var(--text);
      border: 1px solid var(--accent-dim);
      border-radius: 8px;
      cursor: pointer;
      min-width: 180px;
    }}
    select:hover, select:focus {{
      outline: none;
      border-color: var(--accent);
    }}
    .panel {{
      display: none;
      background: var(--card);
      border-radius: 12px;
      padding: 1.5rem;
      box-shadow: 0 4px 6px -1px rgba(0,0,0,0.3);
    }}
    .panel.active {{
      display: block;
    }}
    .family-summary {{
      color: var(--text-muted);
      margin-bottom: 1.5rem;
      font-size: 0.9rem;
    }}
    .member-card {{
      background: var(--bg);
      border-radius: 10px;
      padding: 1.25rem;
      margin-bottom: 1rem;
      border-left: 4px solid var(--accent);
    }}
    .member-card:hover {{
      background: var(--card-hover);
    }}
    .member-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.75rem;
      flex-wrap: wrap;
      gap: 0.5rem;
    }}
    .member-id {{
      font-weight: 600;
      font-size: 1.05rem;
      color: var(--accent);
    }}
    .member-length {{
      color: var(--text-muted);
      font-size: 0.9rem;
    }}
    .bar-container {{
      position: relative;
      height: 24px;
      background: var(--card);
      border-radius: 6px;
      margin-bottom: 0.75rem;
      overflow: hidden;
    }}
    .bar-track {{
      position: absolute;
      left: 0;
      top: 0;
      height: 100%;
      width: 100%;
      background: linear-gradient(90deg, rgba(56,189,248,0.15) 0%, rgba(14,165,233,0.1) 100%);
      border-radius: 6px;
    }}
    .bar-hepn {{
      position: absolute;
      top: 0;
      height: 100%;
      background: var(--hepn);
      opacity: 0.85;
      border-radius: 4px;
    }}
    .hepn-table {{
      width: 100%;
      font-size: 0.85rem;
      border-collapse: collapse;
    }}
    .hepn-table th {{
      text-align: left;
      padding: 0.4rem 0.6rem;
      color: var(--text-muted);
      font-weight: 500;
    }}
    .hepn-table td {{
      padding: 0.4rem 0.6rem;
      border-top: 1px solid var(--card-hover);
    }}
    .no-hepn {{
      color: var(--text-muted);
      font-style: italic;
      font-size: 0.9rem;
    }}
  </style>
</head>
<body>
  <div class="header">
    <h1>Family Dashboard – Deep Hits</h1>
    <p class="subtitle">{len(family_ids)} families · {total_members} total sequences</p>
  </div>

  <div class="controls">
    <label for="family-select">Select family:</label>
    <select id="family-select">
      <option value="">-- Choose family --</option>
      {chr(10).join(f'      <option value="{fid}">{fid}</option>' for fid in family_ids)}
    </select>
  </div>

  <div id="panels"></div>

  <script>
    const data = {json.dumps(data)};

    function renderFamily(fid) {{
      if (!fid || !data[fid]) return '';
      const members = data[fid];
      let html = `<div class="panel active">
        <div class="family-summary">Family ${{fid}}: ${{members.length}} member(s)</div>`;

      members.forEach(m => {{
        html += `<div class="member-card">
          <div class="member-header">
            <span class="member-id">${{m.id}}</span>
            <span class="member-length">${{m.length}} aa</span>
          </div>
          <div class="bar-container">
            <div class="bar-track"></div>`;
        m.hepn.forEach(h => {{
          const left = (h.start / m.length) * 100;
          const w = (h.length / m.length) * 100;
          html += `<div class="bar-hepn" style="left:${{left}}%; width:${{w}}%" title="${{h.motif}}"></div>`;
        }});
        html += `</div>`;
        if (m.hepn.length) {{
          html += `<table class="hepn-table">
            <thead><tr><th>#</th><th>Start</th><th>End</th><th>Length</th><th>Motif</th></tr></thead>
            <tbody>`;
          m.hepn.forEach((h, i) => {{
            html += `<tr><td>${{i+1}}</td><td>${{h.start}}</td><td>${{h.end}}</td><td>${{h.length}}</td><td>${{h.motif}}</td></tr>`;
          }});
          html += `</tbody></table>`;
        }} else {{
          html += `<p class="no-hepn">No HEPN domains detected</p>`;
        }}
        html += `</div>`;
      }});

      html += `</div>`;
      return html;
    }}

    const select = document.getElementById('family-select');
    const panels = document.getElementById('panels');

    function update() {{
      const fid = select.value;
      panels.innerHTML = fid ? renderFamily(fid) : '';
    }}

    select.addEventListener('change', update);
    if (select.value) update();
  </script>
</body>
</html>
'''

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[SUCCESS] Wrote {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate family dashboard from fam_fasta.")
    parser.add_argument(
        "--input",
        default="data/fam_fasta.fasta",
        help="Input FASTA file",
    )
    parser.add_argument(
        "--output",
        default="visualization/family_dashboard.html",
        help="Output HTML file",
    )
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[!] Input file not found: {args.input}")
        return

    families = load_families(args.input)
    if not families:
        print("[!] No families loaded.")
        return

    generate_html(families, args.output)


if __name__ == "__main__":
    main()
