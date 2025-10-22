import itertools as it
from operator import itemgetter
import collections
import plotly.subplots as sp
import plotly.graph_objects as go
import plotly.io as pio
from argparse import ArgumentParser
from PyPDF2 import PdfMerger
import os
from tempfile import TemporaryDirectory

def main():
    args = get_args()
    traces = collections.defaultdict(list)

    for f in args.input:
        sample_full = f.replace(".mosdepth.global.dist.txt", "")
        basename = sample_full.split('/')[11]
        sample = basename.split('_')[0]
        gen = (x.rstrip().split("\t") for x in open(f))
        for chrom, data in it.groupby(gen, itemgetter(0)):
            if chrom.startswith("GL"):
                continue
            if "Un" in chrom: continue
            if "random" in chrom or "HLA" in chrom: continue
            if chrom.endswith("alt"): continue
            xs, ys = [], []
            v50 = 0
            found = False
            for _, x, y in data:
                y = float(y)
                if y < 0.01:
                    continue
                if not found and y > 0.5:
                    v50 = x
                    found = True
                    print(f"{sample}\t{chrom}\t{int(x)}\t{y:.3f}")
                xs.append(float(x))
                ys.append(y)

            if len(xs) > 100:
                xs = [x for i, x in enumerate(xs) if ys[i] > 0.02]
                ys = [y for y in ys if y > 0.02]
                if len(xs) > 100:
                    xs = xs[::2]
                    ys = ys[::2]

            traces[chrom].append({
                'x': xs,
                'y': ys,
                'name': sample + (f" ({float(v50):.1f})" if v50 != 0 else sample)
            })

    chromosomes = list(traces.keys())
    chunk_size = 4
    chunks = [chromosomes[i:i + chunk_size] for i in range(0, len(chromosomes), chunk_size)]

    fig_list = []

    for chunk in chunks:
        fig = sp.make_subplots(
            rows=2, cols=2,
            subplot_titles=[f"Coverage Distribution - {c}" for c in chunk],
            horizontal_spacing=0.15,
            vertical_spacing=0.15
        )
        fig.update_layout(
            height=900,
            width=900,
            margin=dict(l=40, r=40, t=70, b=40),
            showlegend=False,
            autosize=True
        )

        # For each subplot, add data traces and prepare manual legend text/colors
        for i, chrom in enumerate(chunk):
            row = (i // 2) + 1
            col = (i % 2) + 1
            subplot_traces = traces[chrom]
            for trace in subplot_traces:
                fig.add_trace(go.Scatter(
                    x=trace['x'],
                    y=trace['y'],
                    mode='lines',
                    name=trace['name'],
                    line=dict(width=1)
                ), row=row, col=col)
            
            fig.update_xaxes(title_text="Coverage", row=row, col=col)
            fig.update_yaxes(title_text="Proportion of Bases at Coverage", row=row, col=col)

            # Build manual legend text for this subplot
            legend_items = []
            # Pick colors from default Plotly colorway in order
            default_colors = px_default_color_seq()
            for idx, trace in enumerate(subplot_traces):
                color = default_colors[idx % len(default_colors)]
                legend_items.append(f"<span style='color:{color}'>▇</span> {trace['name']}<br>")

            legend_text = "".join(legend_items)

            # Add annotation for this subplot’s legend at top-right corner inside subplot
            # Using domain coordinates relative to subplot (0 to 1)
            x_pos = 0.95  # near right side
            y_pos = 0.95  # near top
            fig.add_annotation(dict(
                x=x_pos,
                y=y_pos,
                xref=f"x{'' if i==0 else i+1} domain",
                yref=f"y{'' if i==0 else i+1} domain",
                text=legend_text,
                showarrow=False,
                align="left",
                font=dict(size=10),
                bordercolor="black",
                borderwidth=1,
                bgcolor="white",
                opacity=0.8,
                xanchor='right',
                yanchor='top'
            ))

        fig_list.append(fig)

    with TemporaryDirectory() as tmpdir:
        pdf_paths = []
        for idx, fig in enumerate(fig_list):
            tmp_pdf = os.path.join(tmpdir, f"page_{idx}.pdf")
            fig.write_image(tmp_pdf, engine="kaleido")
            pdf_paths.append(tmp_pdf)

        merger = PdfMerger()
        for pdf_path in pdf_paths:
            merger.append(pdf_path)

        merger.write(args.output)
        merger.close()

    print(f"Multi-page PDF with separate legends per plot saved to: {args.output}")


def px_default_color_seq():
    # Return the default Plotly qualitative color sequence for line colors
    return [
        "#636efa", "#EF553B", "#00cc96", "#ab63fa", "#FFA15A",
        "#19d3f3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
    ]


def get_args():
    parser = ArgumentParser(description="Creates PDF plots from mosdepth results with 4 plots per page and separate legends per plot.")
    parser.add_argument("-o", "--output",
                        default="dist.pdf",
                        help="Path and name of output PDF file.")
    parser.add_argument("input",
                        nargs='+',
                        help="The .mosdepth.global.dist.txt files to use for plotting.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
