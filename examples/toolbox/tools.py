import matplotlib.pyplot as plt
import numpy as np


def plot_and_compare(mine, *others):
    """Plot and compare the content of a numpy matrix with some others

    - Plot all columns of an input matrix with respect to its first column
    - Compare these columns with the columns of extra args

    e.g.

    >>> ref = np.loadtxt("BouncingBall/BouncingBallTS-MoreauJeanCombinedProjectionOSI.ref", skiprows=1)
    >>> mine = np.loadtxt("BouncingBall/BouncingBallTS-MoreauJeanCombinedProjectionOSI.dat")
    >>> ref2 = np.loadtxt("BouncingBall/BouncingBallTS.ref", skiprows=1)
    >>> plot_and_compare(mine, ref, ref2)

    will plot content of mine, ref and ref2 and (mine-ref), (mine-ref2) ...
    with respect to first column of each matrix.
    """
    num_cols = mine.shape[1]
    num_graphs = len(others) + 1
    colors = ["blue", "green", "orange", "purple", "brown"]
    linestyles = ["-.", ":", "-.", ":"]
    markers = [None, None, None, None, None]
    for i in range(1, num_cols):
        _fig, axes = plt.subplots(
            1, num_graphs, figsize=(6 * num_graphs, 6), sharex=True
        )

        if num_graphs == 1:
            axes = [axes]

        time = mine[:, 0]
        axes[0].plot(
            time, mine[:, i], label=f"Mine - Col {i}", color="red", linestyle="--"
        )

        for j, other in enumerate(others):
            if other is not None and other.shape[1] > i:
                time_other = other[:, 0]
                color = colors[j % len(colors)]
                linestyle = linestyles[j % len(linestyles)]
                marker = markers[j % len(markers)]
                axes[0].plot(
                    time_other,
                    other[:, i],
                    label=f"Other {j + 1} - Col {i}",
                    color=color,
                    linestyle=linestyle,
                    marker=marker,
                )

        axes[0].set_xlabel("Time")
        axes[0].set_title(f"Column {i} vs Time")
        axes[0].legend()
        axes[0].grid()

        # Plot differences for each other matrix
        for j, other in enumerate(others):
            if other is not None and other.shape[1] > i:
                rows = np.min([mine.shape[0], other.shape[0]])
                time_other = other[:rows, 0]
                # y_mine_interp = np.interp(time_other, time, mine[:, i])
                # diff = y_mine_interp - other[:, i]
                diff = mine[:rows, i] - other[:rows, i]
                axes[j + 1].plot(
                    time_other,
                    diff,
                    label=f"Diff Mine - Other {j + 1}",
                    color="red",
                )
                axes[j + 1].set_xlabel("Time")
                axes[j + 1].set_title(f"Diff Col {i} (Mine - Other {j + 1})")
                axes[j + 1].legend()
                axes[j + 1].grid()
