from recognition_data import symmetrySizes, n_distortion_values, max_distortion, shapes, repeats, results, recognizers


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
        line = ", ".join(["{}: {}".format(shapes[i], result_section.count(i))
                          for i in shape_idxs_of_same_size])
        print("{} : {}".format(distortion, line))


if __name__ == "__main__":
    make_table(2, 5)
    make_table(3, 5)
