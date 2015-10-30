#ifndef ERGM_H
#define ERGM_H

typedef struct ergm ergm;

/** Start up embedded R. */
void start_R(void);

/** Stop embedded R. */
void stop_R(void);

/** Initialize an ERGM model.
 *
 * This can later be used to simulate networks with various parameters.
 *
 * \param[in] nnode number of nodes in the model
 * \param[out] an ERGM structure, which can be used for simulations later
 */
ergm *ergm_create(int nnode);

/** Simulate a network from an ERGM.
 *
 * The number of nodes must be the same as that which was passed to
 * ergm_create().
 *
 * \param[out] g uninitialized graph, result will be stored here
 * \param[in] e model to simulate from (see ergm_create())
 * \param[in] nnode number of nodes in the network
 * \param[in] nparam number of ERGM model parameters (for now this must be 2)
 * \param[in] coef vector of parameters to simulate from
 */
void simulate_from_ergm(igraph_t *g, ergm *e, int nnode, int nparam, double *coef);

/** Destroy an ERGM object.
 *
 * \param[in] model ERGM object to destroy
 */
void ergm_free(ergm *model);

#endif
