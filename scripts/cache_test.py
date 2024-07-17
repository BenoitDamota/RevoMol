"""
script to try the cache and joblib

The objective is to use the cache to store results of evaluations for a molecule
instead of recomputing them each time.

To avoid keeping in memory the object Molecule, we have to ignore it in the cache.
Look at the documentation to ignore the object Molecule but keep in mind the SMILES string of the molecule.
https://joblib.readthedocs.io/en/latest/memory.html#ignoring-some-arguments


"""

import time
from functools import cache, wraps
from joblib import Memory


def time_it(func):
    @wraps(func)
    def time_it_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return time_it_wrapper


@time_it
def call_factorial_no_cache(n):
    return factorial_no_cache(n)


@time_it
def call_factorial_with_cache(n):
    return factorial_cache(n)


@time_it
def call_factorial_with_joblib(n):
    return factorial_joblib(n)


def factorial_no_cache(n):
    if n == 0:
        return 1
    time.sleep(1)
    return n * factorial_no_cache(n - 1)


@cache
def factorial_cache(n):
    if n == 0:
        return 1
    time.sleep(1)
    return n * factorial_cache(n - 1)


cachedir = "joblib_cache_test"
memory = Memory(cachedir, verbose=0)


@memory.cache
def factorial_joblib(n):
    if n == 0:
        return 1
    time.sleep(1)
    return n * factorial_joblib(n - 1)


number = 5

print("First call")

call_factorial_no_cache(number)
call_factorial_with_cache(number)
call_factorial_with_joblib(number)

print("Second call")

call_factorial_no_cache(number)
call_factorial_with_cache(number)
call_factorial_with_joblib(number)

print("Clear cache")
factorial_cache.cache_clear()
memory.clear()
print("Third call")


call_factorial_no_cache(number)
call_factorial_with_cache(number)
call_factorial_with_joblib(number)


print("Fourth call")

call_factorial_no_cache(number)
call_factorial_with_cache(number)
call_factorial_with_joblib(number)
