__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from recognition_data import symmetrySizes, n_distortion_values, max_distortion, shapes, repeats, results, recognizers
import matplotlib.pyplot as plt
import numpy as np


def make_table(recognizer_idx, shape_idx):
    print(recognizers[recognizer_idx])
    shape_size = symmetrySizes[shape_idx]
    shape_idxs_of_same_size = [c for c, e in enumerate(
        symmetrySizes) if e == shape_size]
    shape_rec_results = results[(recognizer_idx, shape_idx)]
    for distortion_idx in range(n_distortion_values):
        result_start = distortion_idx * repeats
        result_end = (distortion_idx + 1) * repeats
        result_section = shape_rec_results[result_start:result_end]
        distortion = distortion_idx * \
            max_distortion / (n_distortion_values - 1)
        summed = {shapes[i]: result_section.count(i)
                  for i in shape_idxs_of_same_size}
        line = ", ".join(["{}: {}".format(k, v) for k, v in summed.items()])
        print("{} : {}".format(distortion, line))


def stack_plot(recognizer_idx, shape_idx, axes):
    shape_size = symmetrySizes[shape_idx]
    shape_idxs_of_same_size = [c for c, e in enumerate(
        symmetrySizes) if e == shape_size]
    shape_rec_results = results[(recognizer_idx, shape_idx)]
    sums_per_shape = {k: [] for k in shape_idxs_of_same_size}
    for distortion_idx in range(n_distortion_values):
        result_start = distortion_idx * repeats
        result_end = (distortion_idx + 1) * repeats
        result_section = shape_rec_results[result_start:result_end]
        for i in shape_idxs_of_same_size:
            sums_per_shape[i].append(result_section.count(i))

    cumulative_sums = [0 for i in range(n_distortion_values)]
    assert len(cumulative_sums) == len(sums_per_shape[shape_idx])
    xs = np.arange(n_distortion_values)
    xs_data = np.arange(0.0, 1.1, 0.1)

    # halve number of rows by subsetting everything
    for k, v in sums_per_shape.items():
        sums_per_shape[k] = v[:6]

    cumulative_sums = cumulative_sums[:6]
    xs = xs[:6]
    xs_data = xs_data[:6]

    def autolabel(ax, rects, shape):
        for rect in rects:
            width = rect.get_width()
            if width < 5:
                continue

            annotate_kwargs = {
                "xy": (
                    rect.get_x() + rect.get_width() / 2,
                    rect.get_y() + rect.get_height() / 2
                ),
                "ha": "center",
                "va": "center",
                "color": "white"
            }

            if rect.get_bbox().y0 < 0:
                ax.annotate(shape, weight="bold", **annotate_kwargs)
            else:
                ax.annotate("{}".format(width), **annotate_kwargs)

    for i in shape_idxs_of_same_size:
        rects = axes.barh(y=xs, width=sums_per_shape[i], height=0.8,
                          left=cumulative_sums)
        autolabel(axes, rects, shapes[i])
        cumulative_sums = [a + b for a, b
                           in zip(cumulative_sums, sums_per_shape[i])]

    axes.set_yticks(xs)
    axes.set_yticklabels(["{:.1f}".format(x) for x in xs_data])


if __name__ == "__main__":
    shape_size = 4
    recognizer_labels = ["Ang. & geom. ind.",
                         "CShM", "Biased CShM", "CShM prob. distr."]
    selected_recognizers = [0, 1, 2, 3]
    shape_idxs_of_same_size = [c for c, e in enumerate(symmetrySizes)
                               if e == shape_size]
    fig, axs = plt.subplots(len(selected_recognizers),
                            len(shape_idxs_of_same_size),
                            tight_layout=True, sharey="row")
    for row_idx, recognizer in enumerate(selected_recognizers):
        subplots_row = axs[row_idx]
        subplots_row[0].set_ylabel(recognizer_labels[recognizer])

        for col_idx, shape in enumerate(shape_idxs_of_same_size):
            stack_plot(recognizer, shape, subplots_row[col_idx])

    # Axis hiding
    # Hide y axis for all columns but the first
    for row in axs:
        for plot in row[1:]:
            plot.get_yaxis().set_visible(False)

    # Hide x axis and margins everywhere
    for row in axs:
        for plot in row:
            plot.get_xaxis().set_visible(False)
            plot.margins(0)

    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.show()
