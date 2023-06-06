#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "k.h"

typedef struct geom {
    double sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
} geom; 

void calcgeo(geom *g, double lam, double u1, double u2) {

    (*g).sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
    (*g).cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
    (*g).sigma = atan((*g).sin_sigma / (*g).cos_sigma);
    if ((*g).sigma <= 0) (*g).sigma = M_PI + (*g).sigma;
    (*g).sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / (*g).sin_sigma;
    (*g).cos_sq_alpha = 1 - pow((*g).sin_alpha, 2.);
    (*g).cos2sigma = (*g).cos_sigma - ((2 * sin(u1) * sin(u2)) / (*g).cos_sq_alpha);

};

K vinc(K latpk, K longpk, K latck, K longck) {
    
    geom gvar;

    double req = 6378137.0;             //Radius at equator
    double flat = 1 / 298.257223563;    //flattening of earth
    double rpol = (1 - flat) * req;

    double u1, u2, lon, lam, tol, diff;
    double A, B, C, lam_pre, delta_sig, dis, azi1, usq;

    K res = kf(0);

    if (((latpk->f) + (longpk->f)) == ((latck->f) + (longck->f))) return(res);

    // convert to radians
    double latp = M_PI * (latpk->f) / 180.0;
    double latc = M_PI * (latck->f) / 180.0;
    double longp = M_PI * (longpk->f) / 180.0;
    double longc = M_PI * (longck->f) / 180.0;

    u1 = atan((1 - flat) * tan(latc));
    u2 = atan((1 - flat) * tan(latp));

    lon = longp - longc;
    lam = lon;
    tol = pow(10., -12.); // iteration tolerance
    diff = 1.;

    while (fabs(diff) > tol) {
        calcgeo(&gvar, lam, u1, u2);
        C = (flat / 16) * gvar.cos_sq_alpha * (4 + flat * (4 - 3 * gvar.cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * gvar.sin_alpha * (gvar.sigma + C * gvar.sin_sigma * (gvar.cos2sigma + C * gvar.cos_sigma * (2 * pow(gvar.cos2sigma, 2.) - 1)));
        diff = fabs(lam_pre - lam);
    }
    
    calcgeo(&gvar, lam, u1, u2);
    
    usq = gvar.cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    delta_sig = B * gvar.sin_sigma * (gvar.cos2sigma + 0.25 * B * (gvar.cos_sigma * (-1 + 2 * pow(gvar.cos2sigma, 2.)) -
                                                         (1 / 6) * B * gvar.cos2sigma * (-3 + 4 * pow(gvar.sin_sigma, 2.)) *
                                                         (-3 + 4 * pow(gvar.cos2sigma, 2.))));

    res = kf(0.001 * rpol * A * (gvar.sigma - delta_sig)); // convert to kms
    
    return(res);
}