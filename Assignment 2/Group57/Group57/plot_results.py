#!/usr/bin/env python3
"""
plot_results.py – Generate boxplots of MPI timing results for Assignment 2.

Usage:
    python3 plot_results.py timing_data.csv

Input CSV format (no header):
    <label>,<nx_tag>,<run_id>,<time_seconds>
    e.g.  P32,n120,1,3.456789

The script produces a single PDF file: results_plot.pdf
with two sub-plots (one per data size: nx=120 and nx=240).
"""

import sys
import csv
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_data(path):
    """Return nested dict: data[nx_tag][P_label] = [t1, t2, ...]"""
    data = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(path) as f:
        for row in csv.reader(f):
            row = [r.strip() for r in row]
            if not row or row[0].startswith('#'):
                continue
            p_label, nx_tag, _run, t = row[0], row[1], row[2], float(row[3])
            data[nx_tag][p_label].append(t)
    return data

def make_plot(data, output='results_plot.pdf'):
    p_order   = ['P32', 'P48', 'P64', 'P96']
    p_values  = [32,    48,    64,    96   ]
    nx_tags   = ['n120', 'n240']
    nx_labels = {'n120': 'nx = ny = nz = 120', 'n240': 'nx = ny = nz = 240'}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    fig.suptitle('Assignment 2 – MPI Stencil Timing (d=7, T=5, F=2)', fontsize=13, fontweight='bold')

    for ax, nx_tag in zip(axes, nx_tags):
        series = data.get(nx_tag, {})
        box_data = [series.get(p, [0.0]) for p in p_order]

        bp = ax.boxplot(box_data,
                        patch_artist=True,
                        notch=False,
                        medianprops=dict(color='black', linewidth=2))

        colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_xticks(range(1, len(p_order) + 1))
        ax.set_xticklabels([f'{p}\n(P={v})' for p, v in zip(p_order, p_values)])
        ax.set_xlabel('Number of MPI Processes (P)', fontsize=11)
        ax.set_ylabel('Wall-clock Time (seconds)', fontsize=11)
        ax.set_title(nx_labels[nx_tag], fontsize=11)
        ax.grid(axis='y', linestyle='--', alpha=0.5)
        ax.set_xlim(0.5, len(p_order) + 0.5)

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f'Plot saved to {output}')

if __name__ == '__main__':
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'timing_data.csv'
    data = load_data(csv_path)
    make_plot(data)
