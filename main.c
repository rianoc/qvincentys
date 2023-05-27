#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "k.h"

K vinc(K latpk, K longpk, K latck, K longck) {
    double req = 6378137.0;             //Radius at equator
    double flat = 1 / 298.257223563;    //flattening of earth
    double rpol = (1 - flat) * req;

    double u1, u2, lon, lam, tol, diff, sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
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
        sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
        sigma = atan(sin_sigma / cos_sigma);
        if (sigma <= 0) sigma = M_PI + sigma;
        sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
        cos_sq_alpha = 1 - pow(sin_alpha, 2.);
        cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);
        C = (flat / 16) * cos_sq_alpha * (4 + flat * (4 - 3 * cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * sin_alpha * (sigma + C * sin_sigma * (cos2sigma + C * cos_sigma * (2 * pow(cos2sigma, 2.) - 1)));
        diff = fabs(lam_pre - lam);
    }
    
    sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
    cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
    sigma = atan(sin_sigma / cos_sigma);
    if (sigma <= 0) sigma = M_PI + sigma;
    sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
    cos_sq_alpha = 1 - pow(sin_alpha, 2.);
    cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);

    usq = cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    delta_sig = B * sin_sigma * (cos2sigma + 0.25 * B * (cos_sigma * (-1 + 2 * pow(cos2sigma, 2.)) -
                                                         (1 / 6) * B * cos2sigma * (-3 + 4 * pow(sin_sigma, 2.)) *
                                                         (-3 + 4 * pow(cos2sigma, 2.))));

    res = kf(0.001 * rpol * A * (sigma - delta_sig)); // convert to kms
    
    return(res);
}