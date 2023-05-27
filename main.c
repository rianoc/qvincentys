#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "k.h"

typedef struct geom {
    double sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
} Geom;

Geom geomVar;

struct geom calcGeo(double lam, double u1, double u2) {

    geomVar.sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
    geomVar.cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
    geomVar.sigma = atan(geomVar.sin_sigma / geomVar.cos_sigma);
    if (geomVar.sigma <= 0) geomVar.sigma = M_PI + geomVar.sigma;
    geomVar.sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / geomVar.sin_sigma;
    geomVar.cos_sq_alpha = 1 - pow(geomVar.sin_alpha, 2.);
    geomVar.cos2sigma = geomVar.cos_sigma - ((2 * sin(u1) * sin(u2)) / geomVar.cos_sq_alpha);

    return(geomVar);
};

K vinc(K latpk, K longpk, K latck, K longck) {
    double req = 6378137.0;             //Radius at equator
    double flat = 1 / 298.257223563;    //flattening of earth
    double rpol = (1 - flat) * req;

    double u1, u2, lon, lam, tol, diff;
    double A, B, C, lam_pre, delta_sig, dis, azi1, usq;

    K res = ktn(KF, 2);
    kF(res)[0]=0;
    kF(res)[1]=0;

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
        geomVar = calcGeo(lam, u1, u2);
        C = (flat / 16) * geomVar.cos_sq_alpha * (4 + flat * (4 - 3 * geomVar.cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * geomVar.sin_alpha * (geomVar.sigma + C * geomVar.sin_sigma * (geomVar.cos2sigma + C * geomVar.cos_sigma * (2 * pow(geomVar.cos2sigma, 2.) - 1)));
        diff = fabs(lam_pre - lam);
    }
    
    geomVar = calcGeo(lam, u1, u2);
    
    usq = geomVar.cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    delta_sig = B * geomVar.sin_sigma * (geomVar.cos2sigma + 0.25 * B * (geomVar.cos_sigma * (-1 + 2 * pow(geomVar.cos2sigma, 2.)) -
                                                         (1 / 6) * B * geomVar.cos2sigma * (-3 + 4 * pow(geomVar.sin_sigma, 2.)) *
                                                         (-3 + 4 * pow(geomVar.cos2sigma, 2.))));
    dis = rpol * A * (geomVar.sigma - delta_sig);
    azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));
    
    kF(res)[0]=dis;
    kF(res)[1]=azi1;
    return(res);
}