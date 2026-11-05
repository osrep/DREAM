
#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <gsl/gsl_errno.h>

using namespace std;
using namespace DREAM::FVM;

RadialGridGenerator::RadialGridGenerator(const len_t nr) : nr(nr) {
    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_fmin = gsl_min_fminimizer_alloc(T);
}


RadialGridGenerator::~RadialGridGenerator(){
    gsl_min_fminimizer_free(gsl_fmin);
}


/**
 * Rebuilds magnetic field data and stores all quantities needed for flux surface and bounce averages.
 */
void RadialGridGenerator::RebuildJacobians(RadialGrid *rGrid) {
    real_t realeps = std::numeric_limits<real_t>::epsilon();

    nr = rGrid->GetNr();

    Bmin   = new real_t[GetNr()];
    Bmax   = new real_t[GetNr()];
    Bmin_f = new real_t[GetNr()+1];
    Bmax_f = new real_t[GetNr()+1];
    theta_Bmin     = new real_t[GetNr()];
    theta_Bmax     = new real_t[GetNr()];
    theta_Bmin_f   = new real_t[GetNr()+1];
    theta_Bmax_f   = new real_t[GetNr()+1];
    xi0TrappedBoundary   = new real_t[GetNr()];
    xi0TrappedBoundary_f = new real_t[GetNr()+1];

    for (len_t ir = 0; ir < GetNr(); ir++){
        theta_Bmin[ir] = getTheta_Bmin(ir);
        theta_Bmax[ir] = getTheta_Bmax(ir);
        if(isUpDownSymmetric && fabs(fabs(theta_Bmax[ir]-theta_Bmin[ir])-M_PI) > 1e-3)
            throw FVMException("RadialGridGenerator: theta_B_max - theta_B_min != pi in an up-down-symmetric B field."
                " Either shaping (delta) and/or poloidal fields (psiPrimeRef) are too strong.");

        Bmin[ir] = BAtTheta(ir,theta_Bmin[ir]);
        Bmax[ir] = BAtTheta(ir,theta_Bmax[ir]);
        if(!Bmin[ir] || 1-Bmin[ir]/Bmax[ir]<100*realeps)
            xi0TrappedBoundary[ir] = 0;
        else
            xi0TrappedBoundary[ir] = sqrt(1-Bmin[ir]/Bmax[ir]);

    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        theta_Bmin_f[ir] = getTheta_Bmin_f(ir);
        theta_Bmax_f[ir] = getTheta_Bmax_f(ir);
        if(isUpDownSymmetric && fabs(fabs(theta_Bmax_f[ir]-theta_Bmin_f[ir])-M_PI) > 1e-3)
            throw FVMException("RadialGridGenerator: theta_B_max - theta_B_min != pi in an up-down-symmetric B field."
                " Either shaping (delta) and/or poloidal fields (psiPrimeRef) are too strong.");
        Bmin_f[ir] = BAtTheta_f(ir,theta_Bmin_f[ir]);
        Bmax_f[ir] = BAtTheta_f(ir,theta_Bmax_f[ir]);
        if(!Bmin_f[ir] || 1-Bmin_f[ir]/Bmax_f[ir]<100*realeps)
            xi0TrappedBoundary_f[ir] = 0;
        else
            xi0TrappedBoundary_f[ir] = sqrt(1-Bmin_f[ir]/Bmax_f[ir]);
    }
    rGrid->SetMagneticExtremumData(
        Bmin, Bmin_f, Bmax, Bmax_f, 
        theta_Bmin, theta_Bmin_f, theta_Bmax, theta_Bmax_f,
        xi0TrappedBoundary, xi0TrappedBoundary_f
    );
}


// Evaluates the magnetic field strength at radial index ir 
// on the distribution grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta(const len_t ir, const real_t theta) {
    real_t ROverR0 = ROverR0AtTheta(ir,theta);
    real_t Btor = BtorGOverR0[ir]/ROverR0;
    real_t Bpol = 0;
    if(psiPrimeRef[ir])
        Bpol = sqrt(NablaR2AtTheta(ir,theta))*psiPrimeRef[ir]/(2*M_PI*ROverR0);  
    return sqrt(Btor*Btor+Bpol*Bpol);
}
// Evaluates the magnetic field strength at radial index ir 
// on the radial flux grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta_f(const len_t ir, const real_t theta) {
    real_t ROverR0 = ROverR0AtTheta_f(ir,theta);
    real_t Btor = BtorGOverR0_f[ir]/ROverR0;
    real_t Bpol = 0;
    if(psiPrimeRef_f[ir])
        Bpol = sqrt(NablaR2AtTheta_f(ir,theta))*psiPrimeRef_f[ir]/(2*M_PI*ROverR0);  
    return sqrt(Btor*Btor+Bpol*Bpol);
}

// The remaining functions are related to determining theta_Bmin and theta_Bmax
// with a gsl fmin algorithm
struct EvalBParams {len_t ir; RadialGridGenerator* rgg; int_t sgn;};
real_t gslEvalB(real_t theta, void *par){
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtTheta(params->ir,theta);
}
real_t gslEvalB_f(real_t theta, void *par){
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtTheta_f(params->ir,theta);
}


// Tolerances and max number of iterations for gsl magnetic field minimizer
const len_t MAX_NUM_ITER = 30;
const real_t EPSABS = 1e-6;
const real_t EPSREL = 0;
/**
 * Finds the extremum of the magnetic field on the interval [0,pi]. 
 * If sgn=1, returns the minimum.
 * If sgn=-1, returns the maximum.
 */
real_t RadialGridGenerator::FindMagneticFieldExtremum(
    len_t ir, int_t sgn, fluxGridType fluxGridType
) {
    real_t theta = FindMagneticFieldExtremumInHalfPlane(ir, sgn, fluxGridType);

    // The search above covers one half plane, [0, pi] or [pi, 2*pi]. When it
    // returns an edge of that half plane (theta = 0 or pi), the extremum can lie
    // just across the edge, in the other half plane: the field then keeps
    // decreasing (for sgn = 1) on one side of the edge. On the ITER LUKE
    // equilibrium the minimum of B lies at theta = 0.0216 on some surfaces while
    // the search returned theta = 0, with B(0.0216)/B(0) - 1 = -7.3e-6. Check
    // both sides of an edge and search the half plane that continues downhill.
    real_t edge;
    if (theta < 2*EPSABS || theta > 2*M_PI - 2*EPSABS)
        edge = 0;
    else if (fabs(theta - M_PI) < 2*EPSABS)
        edge = M_PI;
    else
        return theta;

    EvalBParams params = {ir, this, sgn};
    gsl_function gsl_func;
    gsl_func.function = (fluxGridType == FLUXGRIDTYPE_RADIAL) ? &(gslEvalB_f) : &(gslEvalB);
    gsl_func.params = &(params);

    const real_t d = 10*EPSABS;
    const real_t f_edge = gsl_func.function(edge, gsl_func.params);
    const real_t f_above = gsl_func.function(edge + d, gsl_func.params);
    const real_t f_below = gsl_func.function(edge == 0 ? 2*M_PI - d : edge - d, gsl_func.params);
    if (f_above < f_edge && f_above <= f_below)
        return MinimizeOnInterval(&gsl_func, edge + d, edge, edge + M_PI);
    else if (f_below < f_edge) {
        if (edge == 0)
            return MinimizeOnInterval(&gsl_func, 2*M_PI - d, M_PI, 2*M_PI);
        else
            return MinimizeOnInterval(&gsl_func, M_PI - d, 0, M_PI);
    }
    return edge;
}

/**
 * Minimizes gsl_func on [theta_lower, theta_upper], starting from theta_guess,
 * which must have a lower function value than both interval limits. Returns a
 * limit when the minimum lies within 2*EPSABS of it, and 0 for 2*pi.
 */
real_t RadialGridGenerator::MinimizeOnInterval(
    gsl_function *gsl_func, real_t theta_guess, real_t theta_lower, real_t theta_upper
) {
    // Without a bracket the minimizer cannot start: keep the limit next to the guess.
    const real_t f_guess = gsl_func->function(theta_guess, gsl_func->params);
    if (!(f_guess < gsl_func->function(theta_lower, gsl_func->params))
        || !(f_guess < gsl_func->function(theta_upper, gsl_func->params))) {
        real_t limit = (theta_guess - theta_lower < theta_upper - theta_guess) ? theta_lower : theta_upper;
        return (limit >= 2*M_PI) ? limit - 2*M_PI : limit;
    }
    gsl_min_fminimizer_set(gsl_fmin, gsl_func, theta_guess, theta_lower, theta_upper);
    for (len_t iter=0; iter<MAX_NUM_ITER; iter++) {
        gsl_min_fminimizer_iterate(gsl_fmin);
        real_t x_lo = gsl_min_fminimizer_x_lower(gsl_fmin);
        real_t x_up = gsl_min_fminimizer_x_upper(gsl_fmin);
        if (gsl_min_test_interval(x_lo, x_up, EPSABS, EPSREL) == GSL_SUCCESS)
            break;
    }
    real_t extremum = gsl_min_fminimizer_x_minimum(gsl_fmin);
    if (fabs(extremum - theta_lower) < 2*EPSABS)
        extremum = theta_lower;
    else if (fabs(extremum - theta_upper) < 2*EPSABS)
        extremum = theta_upper;
    return (extremum >= 2*M_PI) ? extremum - 2*M_PI : extremum;
}

/**
 * Finds the extremum of the magnetic field in one half plane, [0, pi] or
 * [pi, 2*pi], chosen from a guess. See FindMagneticFieldExtremum().
 */
real_t RadialGridGenerator::FindMagneticFieldExtremumInHalfPlane(
    len_t ir, int_t sgn, fluxGridType fluxGridType
) {
    real_t theta_lim_lower = 0, theta_lim_upper = M_PI;
	real_t theta_guess;

	EvalBParams params = {ir, this, sgn};
	gsl_function gsl_func;
	if(fluxGridType == FLUXGRIDTYPE_RADIAL)
		gsl_func.function = &(gslEvalB_f);
	else
		gsl_func.function = &(gslEvalB);
	gsl_func.params = &(params);

	if (this->CanGuessThetaOptimum()) {
		if (sgn == 1)
			theta_guess = this->GetThetaBminGuess(ir, fluxGridType);
		else
			theta_guess = this->GetThetaBmaxGuess(ir, fluxGridType);

		/*
		// Adjust guess to lie within correct intervall
		if (theta_guess < theta_lim_lower)
			theta_guess += 2*M_PI;
		else if(theta_guess > theta_lim_upper)
			theta_guess -= 2*M_PI;
		*/

		// Adjust search interval so that 'theta_guess' lies within it
		if (theta_guess > M_PI)
			theta_lim_lower = M_PI, theta_lim_upper = 2*M_PI;

		// If the guess coincides with an interval limit, an
		// error will be thrown by GSL...
		if (theta_guess == theta_lim_lower)
			theta_guess = theta_lim_lower + 10*EPSABS;
		else if (theta_guess == theta_lim_upper)
			theta_guess = theta_lim_upper - 10*EPSABS;
	} else {
		// Adjust the limits if the magnetic field is not up-down symmetric...
		if(!isUpDownSymmetric) {
			real_t
				B, Beps, Bg, guess = sgn>0 ? 0 : M_PI;

			// Evaluate magnetic field on both sides of the guess point
			// to determine in which direction (upper/lower half plane)
			// the extremum lies...
			if (fluxGridType == FLUXGRIDTYPE_DISTRIBUTION) {
				Bg = BAtTheta(ir, guess);
				B = BAtTheta(ir, guess+sgn*1e-3);
				Beps = BAtTheta(ir, 2*M_PI-(guess+sgn*1e-3));
			} else {
				Bg = BAtTheta_f(ir, guess);
				B = BAtTheta_f(ir, guess+sgn*1e-3);
				Beps = BAtTheta_f(ir, 2*M_PI-(guess+sgn*1e-3));
			}

			// Is extremum in upper half plane?
			if (((sgn*Bg < sgn*B) && (sgn*Bg < sgn*Beps)) || sgn*B < sgn*Beps)
				theta_lim_lower = 0, theta_lim_upper = M_PI;
			else
				theta_lim_lower = M_PI, theta_lim_upper = 2*M_PI;
		}
	}

	if((sgn==1 && theta_lim_lower < M_PI) || (sgn==-1 && theta_lim_upper > M_PI)) {
		// if B has a local minimum in theta=theta_lower, return theta_lower
		theta_guess = theta_lim_lower + 10*EPSABS;
		if(gsl_func.function(theta_guess,gsl_func.params) >= gsl_func.function(theta_lim_lower,gsl_func.params))
			return theta_lim_lower;
	} else {
		// if B has a local maximum in theta=pi, return pi
		theta_guess=theta_lim_upper-10*EPSABS;
		if(gsl_func.function(theta_guess,gsl_func.params) >= gsl_func.function(theta_lim_upper,gsl_func.params)) {
			if (theta_lim_upper == 2*M_PI)
				return 0;
			else
				return theta_lim_upper;
		}
	}

    // find extremum with fmin algorithm
    gsl_min_fminimizer_set(
        gsl_fmin, &gsl_func, theta_guess, theta_lim_lower, theta_lim_upper
    );

    int status;
    for(len_t iter=0; iter<MAX_NUM_ITER; iter++){
        gsl_min_fminimizer_iterate (gsl_fmin);
        real_t x_lo = gsl_min_fminimizer_x_lower (gsl_fmin);
        real_t x_up = gsl_min_fminimizer_x_upper (gsl_fmin);
        status = gsl_min_test_interval(x_lo, x_up, EPSABS, EPSREL);
        if(status == GSL_SUCCESS)
            break;
    }
    real_t extremum = gsl_min_fminimizer_x_minimum(gsl_fmin); 
    if(extremum < 2*EPSABS)
        return 0;
    else if (fabs(M_PI-extremum) < 2*EPSABS)
        return M_PI;
    else
        return extremum; 
}

/**
 * Finds the minor radius coordinate of the point of closest approach to the magnetic axis 
 * along the line between (x1,y1,z1) and (x2,y2,z2)
 *
 * Here, this is done by solving dr/dtau = 0 using a bisection scheme, where 0<tau<1 
 * parameterises the line between (x1,y1,z1) and (x2,y2,z2).
 */
real_t RadialGridGenerator::FindClosestApproach(
    real_t x1, real_t y1, real_t z1,
    real_t x2, real_t y2, real_t z2
) {

    // We use the total distance between (x1,y1,z1) and (x2,y2,z2) as a characteristic
    // length scale in this context
    real_t lengthScale=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    
    real_t r;
    real_t theta;
    real_t phi;
    real_t *gradr = new real_t[3];
    real_t drdtau;
    real_t x_tmp;
    real_t y_tmp;
    real_t z_tmp;
    do{
        x_tmp=(x1+x2)/2;
        y_tmp=(y1+y2)/2;
        z_tmp=(z1+z2)/2;
        this->GetRThetaPhiFromCartesian(&r, &theta, &phi, x_tmp, y_tmp, z_tmp, lengthScale, r);
        this->GetGradRCartesian(gradr, r, theta, phi);
        drdtau = gradr[0]*(x2-x1) + gradr[1]*(y2-y1) + gradr[2]*(z2-z1);
        if(drdtau<0){
            x1 = x_tmp;
            y1 = y_tmp;
            z1 = z_tmp;
        }else{
            x2 = x_tmp;
            y2 = y_tmp;
            z2 = z_tmp;
        }
    }while(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))>lengthScale*CartesianCoordinateTol);
    
    delete [] gradr;
    
    return r;
}

