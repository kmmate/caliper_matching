# Caliper Matching Estimators of Average Treatment Effects

C-implementation of the caliper matching estimator. See: arXiv paper [Kormos, Van der Pas, Van der Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects](https://arxiv.org/abs/2304.08373).

## Organisation

- `cm.h` and `cm.c`: main functionalities of the library, a user interface for the caliper matching implementation.
- `example_simple.c`: illustrates the use of the library with two simple examples.
- `example_simulation.c`: illustrates the use of the library in a Monte Carlo simulation setting.

The remaining files in `src` are simply supporting ones.

## Workflow Outline

Illustration of the main steps in the workflow. For a working example, refer to `src/example_simple.c`.

0. Include headers:

        #include <stdio.h>
        #include <stdlib.h>
        #include "cm.h"

1. Choose between models corresponding to either known or estimated propensity score. For estimated propensity score, allocate the model with
        
        CMModel *cm_model = malloc(sizeof(CMModel));

2. Pass data and set estimation settings (assuming `y`, `d`, `x`, `theta` are defined).

        cm_model->y = y; // outcome
        cm_model->d = d; // treatment
        cm_model->x = x; // covariates
        cm_model->delta = 0.0; // use default data-driven caliper
        cm_model->modeltype = "logit";  // propensity score model
        cm_model->theta = theta;  // estimated propensity score parameter
        // ... other settings

3. Initialise the model.
        
        cm_initialise(cm_model); // performs tests on inputs, computes caliper

4. Perform estimation.

        CMResults *results = malloc(sizeof(CMResults)); // to store estimation results
        results = cm_cm(cm_model);  // estimation

5. View results.

        printf("ate_hat = %f.\n", results->ate_hat);
        printf("att_hat = %f.\n", results->att_hat);




## Notes

- The implementation depends on [GSL](https://www.gnu.org/software/gsl/) for linear algebra and on `pthreads` for parallel execution.

- On a Mac, `example_simple.c` can be compiled by
        
        $ gcc -pthread -lgsl -l gslcblas example.c cm.c vector.c matrix.c dynamicarray.c linkedlist.c propscore.c nonpara.c linalg.c -o example_simple

and then run the binary as
        
        $ ./example_simple
.

- For appropriate parallelisation, adjust the `NUM_THREADS` definition in `cm.h`.