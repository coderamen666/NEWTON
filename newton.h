/* NEWTON
	A RK4 solver for the integration of Newton's equations

	Written by Jackie R.
	Use this however you want (preferrably not for evil!)
*/

// NOTE: This header file, newton.h, is only useful
// in conjunction with the implementation in newton.c

#ifndef _NEWTON_H_
#define _NEWTON_H_
#include <stdlib.h>

// HERE LIES THE BEGINNING OF THE DOCUMENTATION AND FUNCTION/STRUCT DECLARATIONS
// DON'T CHANGE THESE UNLESS YOU KNOW WHAT YOU'RE DOING!

/*
	newton_solution_t - data structure that represents the output of the RK4 solver

	Elements:
		len - number of indexed sample points of the solution
		t - array of times that correspond to integer index
		v - array of derivatives that correspond to indexed times
		x - array of function values that correspond to indexed times
*/

typedef struct {
	int len;
	double *t;
	double *v;
	double *x;
} newton_solution_t;

/*
	newton_solve - solves non-coupled 2nd order differential equations using RK4
		(like those found in Newtonian mechanics!)

	Arguments:
		a - function for the acceleration
			Its signature is double a(double x, double v, double t);
		x0, v0, t0 - Initial position, velocity, and time, respectively
		dt - the length of a given timestep
		steps - the number of timesteps to take during the solution

	Returns:
		A pointer to a freshly allocated newton_solution_t (or NULL if the allocation fails)
*/
newton_solution_t *newton_solve(double (*a)(double,double,double), \
	double x0, double v0, double t0, double dt, unsigned int steps);

/*
	newton_dealloc_solution - function to deallocate pointers to newton_solution_t

	Arguments:
		sol - pointer to newton_solution_t

	Returns:
		nothing
*/
void newton_dealloc_solution(newton_solution_t *sol);

#endif
