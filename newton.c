/* NEWTON
	A RK4 solver for the integration of Newton's equations

	Written by Jackie R.
	Use this however you want (preferrably not for evil!)
*/

#include "newton.h"

// The big, scary solver. This function is the meat of the package
newton_solution_t *newton_solve(double (*a)(double,double,double), \
	double x0, double v0, double t0, double dt, unsigned int steps){

	// Allocate memory
	newton_solution_t *sol = calloc(1, sizeof(newton_solution_t));

	double *v = calloc(steps + 1, sizeof(double));
	double *x = calloc(steps + 1, sizeof(double));
	double *t = calloc(steps + 1, sizeof(double));

	if(!sol || !v || !x || !t) {
		// The allocation has failed!
		return NULL;
	}

	// Temporary variables
	double kv[4];
	double kx[4];

	// Set initial values
	x[0] = x0;
	v[0] = v0;

	// Solver loop
	for(unsigned int i = 0; i < steps; i++) {
		// Update time
		t[i] = (dt * i) + t0;

		// Calculate intermediary slopes
		kx[0] = v[i] * dt;
		kv[0] = (*a)(x[i], v[i], t[i]) * dt;
		kx[1] = (v[i] + .5*kv[0]) * dt;
		kv[1] = (*a)(x[i] + .5*kx[0], v[i] + .5*kv[0], t[i] + .5*dt) * dt;
		kx[2] = (v[i] + .5*kv[1]) * dt;
		kv[2] = (*a)(x[i] + .5*kx[1], v[i] + .5*kv[1], t[i] + .5*dt) * dt;
		kx[3] = (v[i] + kv[2]) * dt;
		kv[3] = (*a)(x[i] + kx[2], v[i] + kv[2], t[i] + dt) * dt;

		// Update dynamical variables
		x[i + 1] = x[i] + (1.0/6)*(kx[0] + 2*kx[1] + 2*kx[2] + kx[3]);
		v[i + 1] = v[i] + (1.0/6)*(kv[0] + 2*kv[1] + 2*kv[2] + kv[3]);
	}
	// Update that last time
	t[steps] = (dt * steps) + t0;

	// Pack up the newton_solution_t and return the approximate solution
	sol->t = t;
	sol->x = x;
	sol->v = v;
	sol->len = steps + 1;

	return sol;
}

void newton_dealloc_solution(newton_solution_t *sol) {
	// Liberating your computer's RAM from the tyranny of Issac Newton!
	free(sol->x);
	free(sol->v);
	free(sol->t);
	free(sol);
}
