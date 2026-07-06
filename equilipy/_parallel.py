"""Shared process-based parallel execution helpers."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import Any, TypeVar

from joblib import Parallel, delayed
from tqdm import tqdm

T = TypeVar("T")


def starmap_joblib(
    function: Callable[..., T],
    args: Iterable[Sequence[Any]],
    n_jobs: int | None,
    *,
    progress: bool = False,
    progress_callback: Callable[[int, int, str], None] | None = None,
    **progress_kwargs: Any,
) -> list[T]:
    """Run a starmap-style workload with joblib processes."""
    arg_list = list(args)
    if not arg_list:
        return []

    worker_count = 1 if n_jobs is None else int(n_jobs)
    if worker_count < 1:
        worker_count = 1
    total = len(arg_list)
    label = str(progress_kwargs.get("desc") or "Batch")

    def notify(completed: int) -> None:
        if progress_callback is not None:
            progress_callback(completed, total, label)

    if progress and progress_callback is None:
        arg_iterator = tqdm(arg_list, **progress_kwargs)
    else:
        arg_iterator = arg_list

    if worker_count == 1:
        notify(0)
        results = []
        for index, arg in enumerate(arg_iterator, start=1):
            results.append(function(*arg))
            notify(index)
        return results

    if progress_callback is None:
        with Parallel(
            n_jobs=worker_count,
            backend="loky",
            max_nbytes=None,
        ) as parallel:
            return parallel(delayed(function)(*arg) for arg in arg_iterator)

    notify(0)
    results = []
    with Parallel(
        n_jobs=worker_count,
        backend="loky",
        max_nbytes=None,
        return_as="generator",
    ) as parallel:
        result_iterator = parallel(delayed(function)(*arg) for arg in arg_list)
        for index, result in enumerate(result_iterator, start=1):
            results.append(result)
            notify(index)
    return results
