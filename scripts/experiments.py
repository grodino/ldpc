from pathlib import Path

import json
from typing import IO, Optional
import click
import numpy as np
from scipy.special import erfc
from matplotlib import pyplot as pl

# Setup the figures directory
SCRIPT_DIR = Path(__file__).resolve().parent

png_fig = SCRIPT_DIR.joinpath('figures/')
pgf_fig = SCRIPT_DIR.parent.joinpath('report/figures/')
png_fig.mkdir(exist_ok=True, parents=True)
pgf_fig.mkdir(exist_ok=True, parents=True)

TEXT_WIDTH = 5.78851  # inches
pl.rcParams.update({
    "figure.figsize": [.6*6.4, .6*4.8],     # change figure default size
    "savefig.bbox": "tight",                # image fitted to the figure
    # grid lines for major and minor ticks
    "axes.grid.which": "major",
    "axes.grid": True,                      # always display the grid
    "lines.marker": "+",                    # add a marker to each point
    "lines.linewidth": .7,                  # reduce linewidth to better see the points
    "font.family": "serif",                 # use serif/main font for text elements
    "text.usetex": True,                    # use inline math for ticks
    "pgf.rcfonts": False,                   # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        r"\usepackage{url}",                # load additional packages
        r"\usepackage{unicode-math}",       # unicode math setup
        r"\setmainfont{DejaVu Serif}",      # serif font via preamble
    ])
})


@click.command()
@click.argument('input', type=click.File('r'))
@click.option('--save', type=bool, is_flag=True)
def main(input: IO, save: bool):
    experiments = input.read(-1)

    name = Path(input.name).stem
    if not(experiments):
        raise ValueError('Did not receive information')

    experiments: dict = json.loads(experiments)

    for codes in experiments['codes']:
        figure, (bler_ax, ber_ax) = pl.subplots(1, 2)
        figure.set_size_inches(TEXT_WIDTH, TEXT_WIDTH / 2.3)
        figure.tight_layout(w_pad=2)

        for code in codes:
            noise_parameter = np.array(code['noise_parameter'])

            ber_ax.semilogy(
                10 * np.log10(noise_parameter),
                code['bit_error_rate'],
                '-+',
                label=f"${code['name']}$",
                nonpositive='mask'
            )
            bler_ax.semilogy(
                10 * np.log10(noise_parameter),
                code['block_error_rate'],
                '-+',
                label=f"${code['name']}$",
                nonpositive='mask'
            )

        uncoded_ber = erfc(
            np.sqrt(2 * noise_parameter)
        )
        ber_ax.semilogy(
            10 * np.log10(noise_parameter),
            uncoded_ber,
            '-+',
            label='uncoded'
        )

        code_len = len(code['code']['parity_col'])
        uncoded_bler = 1 - (1 - uncoded_ber)**code_len
        bler_ax.semilogy(
            10 * np.log10(noise_parameter),
            uncoded_bler,
            '-+',
            label='uncoded'
        )

        bler_ax.axvline(10 * np.log10(
            (2**code['code_rate'] - 1) /
            (2 * code['code_rate'])
        ))
        ber_ax.axvline(10 * np.log10(
            (2**code['code_rate'] - 1) /
            (2 * code['code_rate'])
        ))

        bler_ax.set_xlabel('SNR (dB)')
        bler_ax.set_ylabel('$P_b$')
        ber_ax.set_xlabel('SNR (dB)')
        ber_ax.set_ylabel('$P_e$')

        handles, labels = ber_ax.get_legend_handles_labels()
        figure.legend(
            handles,
            labels,
            bbox_to_anchor=(-.03, .95, 1, 0),
            loc='lower left',
            ncol=5,
            columnspacing=.5,
            mode='expand',
            handlelength=1,
        )

        if save:
            figure.savefig(png_fig.joinpath(
                f'{name}_rate-{code["code_rate"]:0.2}.png'
            ), dpi=600)
            figure.savefig(pgf_fig.joinpath(
                f'{name}_rate-{code["code_rate"]:0.2}.pgf'
            ), dpi=600)

    if not(save):
        pl.legend()
        pl.show()


if __name__ == '__main__':
    main()
