from __future__ import annotations

from joblib import Parallel, delayed
from tqdm import tqdm

import warnings


def configure_worker_warnings():
    """
    Configure warnings to ignore in parallel workers.
    """

    warnings.filterwarnings(
        "ignore",
        message=r"CrystalNN: cannot locate an appropriate radius.*",
        category=UserWarning,
    )

    warnings.filterwarnings(
        "ignore",
        message=r"No oxidation states specified on sites!*",
        category=UserWarning,
    )

    warnings.filterwarnings(
        "ignore",
        message=r"Issues encountered while parsing CIF:*",
        category=UserWarning,
    )


def parallel_map(
    items,
    function,
    workers,
    description,
):
    """
    Apply a function to items in parallel with a progress bar.

    Parameters
    ----------
    items
        Iterable of input items.

    function
        Function executed for each item.

    workers
        Number of parallel workers.

    description
        Progress-bar description.

    Returns
    -------
    list
        Results returned by the worker function.
    """

    def worker(item):
        configure_worker_warnings()
        return function(item)

    iterator = tqdm(
        items,
        desc=description,
    )

    return Parallel(
        n_jobs=workers,
    )(
        delayed(worker)(item)
        for item in iterator
    )
