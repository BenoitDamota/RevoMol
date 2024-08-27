"""
script to try the cache and joblib

The objective is to use the cache to store results of evaluations for a molecule
instead of recomputing them each time.

To avoid keeping in memory the object Molecule, we have to ignore it in the cache.
Look at the documentation to ignore the object Molecule but keep in mind the 
SMILES string of the molecule.
https://joblib.readthedocs.io/en/latest/memory.html#ignoring-some-arguments


"""

import time
from functools import cache, wraps
from typing import Callable

from joblib import Memory


def time_it(func: Callable[[int], int]) -> Callable[[int], int]:
    """Decorator to time the execution of a function"""

    @wraps(func)
    def time_it_wrapper(*args, **kwargs):  # type: ignore[no-untyped-def]
        """Wrapper to time the execution of a function"""
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return time_it_wrapper


@time_it
def call_factorial_no_cache(n: int) -> int:
    """call the factorial function without cache"""
    return factorial_no_cache(n)


@time_it
def call_factorial_with_cache(n: int) -> int:
    """call the factorial function with cache"""
    return factorial_cache(n)


@time_it
def call_factorial_with_joblib(n: int) -> int:
    """call the factorial function with joblib cache"""
    return int(factorial_joblib(n))


def factorial_no_cache(n: int) -> int:
    """factorial function without cache"""
    if n == 0:
        return 1
    time.sleep(1)
    return n * factorial_no_cache(n - 1)


@cache
def factorial_cache(n: int) -> int:
    """factorial function with cache"""
    if n == 0:
        return 1
    time.sleep(1)
    return n * factorial_cache(n - 1)


CACHE_DIR = "joblib_cache_test"
memory = Memory(CACHE_DIR, verbose=0)


@memory.cache  # type: ignore[misc]
def factorial_joblib(n: int) -> int:
    """factorial function with joblib cache"""
    if n == 0:
        return 1
    time.sleep(1)
    return n * int(factorial_joblib(n - 1))


NUMBER = 5

print("First call")

call_factorial_no_cache(NUMBER)
call_factorial_with_cache(NUMBER)
call_factorial_with_joblib(NUMBER)

print("Second call")

call_factorial_no_cache(NUMBER)
call_factorial_with_cache(NUMBER)
call_factorial_with_joblib(NUMBER)

print("Clear cache")
factorial_cache.cache_clear()
memory.clear()
print("Third call")


call_factorial_no_cache(NUMBER)
call_factorial_with_cache(NUMBER)
call_factorial_with_joblib(NUMBER)


print("Fourth call")

call_factorial_no_cache(NUMBER)
call_factorial_with_cache(NUMBER)
call_factorial_with_joblib(NUMBER)
