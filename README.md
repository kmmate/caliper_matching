# Caliper Matching Estimators of Average Treatment Effects

C-implementation of the caliper matching estimator. See: arXiv paper [Kormos, Van der Pas, Van der Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects](https://arxiv.org/abs/2304.08373).

## Organisation

- `cm.h` and `cm.c`: main functionalities of the library, a user interface for the caliper matching implementation.
- `example_simple.c`: illustrates the use of the library with two simple examples.
- `example_simulation.c`: illustrates the use of the library in a Monte Carlo simulation setting.

The remaining files in `src` are simply supporting ones.

## Notes

- The implementation depends on [GSL](https://www.gnu.org/software/gsl/) for linear algebra and on `pthreads` for parallel execution.

- On a Mac, `example_simple.c` can be compiled by
    gcc -pthread -lgsl -l gslcblas example.c cm.c vector.c matrix.c dynamicarray.c linkedlist.c propscore.c nonpara.c linalg.c -o example_simple
and then run the binary as
    ./example_simple
.

- For appropriate parallelisation, adjust the `NUM_THREADS` definition in `cm.h`.