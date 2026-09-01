import pandas as pd
import sys 

species_file = sys.argv[1]
color_file = sys.argv[2]
output_html = sys.argv[3]

species_df = pd.read_csv(species_file, sep='\t', header=None,
                         names=["domain", "phylum", "order", "family", "species"])
with open(color_file) as f:
    colors = [line.strip() for line in f if line.strip()]

species_df["color"] = colors[:len(species_df)]

html = [
    "<!DOCTYPE html>",
    "<html><head><title>Color Preview with Species</title>",
    "<style>body { font-family: monospace; font-size: 12px; } .color-row { display: flex; align-items: center; margin: 2px 0; } .color-box { width: 20px; height: 20px; min-width: 20px; margin-right: 10px; border: 1px solid #000; display: inline-block; } h2 { margin-top: 20px; border-bottom: 2px solid #333; } h3 { margin: 15px 0 5px 0; color: #555; }</style>",
    "</head><body>",
    "<h1>Species Color Preview</h1>"
]

current_domain = None
current_order = None

for _, row in species_df.iterrows():
    if row["domain"] != current_domain:
        html.append(f"<h2>{row['domain']}</h2>")
        current_domain = row["domain"]
        current_order = None
    if row["order"] != current_order:
        html.append(f"<h3>{row['order']}</h3>")
        current_order = row["order"]
    # Use the exact species name from the last column
    html.append(
        f'<div class="color-row"><div class="color-box" style="background-color:#{row["color"]};"></div><span>{row["species"]} - #{row["color"]}</span></div>'
    )

html.append("</body></html>")

with open(output_html, "w") as f:
    f.write("\n".join(html))

print(f"Updated {output_html} based on new species order and colors, using exact species names from the file.")