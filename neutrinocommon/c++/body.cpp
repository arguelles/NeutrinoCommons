#include "body.h"
#include "physconst.h"
#include <math.h>

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

PhysConst param;

/*
-----------------------------------------------------------------------
         BODY CLASS DEFINITIONS
-----------------------------------------------------------------------         
*/

Body::Body()
        {
            name = "GenericBody";
        }        
        
// track constructor        
Body::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Body::Track::Track()
        {
            //pass
        }        

double Body::density(Track track_input){
            return 0.0;
        }
        
double Body::ye(Track track_input){
            return 1.0;
        }
        
/*
----------------------------------------------------------------------
         VACUUM CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// vacuum constructor
Vacuum::Vacuum(string name_input = "Vacuum")
        {
            name = name_input;
        }
Vacuum::Vacuum()
        {
            name = "Vacuum";
        }        
        
// track constructor        
Vacuum::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Vacuum::Track::Track()
        {
            //pass
        }        

double Vacuum::density(Track track_input){
            return 0.0;
        }
        
double Vacuum::ye(Track track_input){
            return 1.0;
        }
        
/*
----------------------------------------------------------------------
         ConstantDensity CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
ConstantDensity:ConstantDensity(density_input,ye_input)
        {
            name = "ConstantDensity";
            constant_density = density_input;
            constant_ye = ye_input;
        }
// track constructor        
ConstantDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
ConstantDensity::Track::Track()
        {
            //pass
        }
double ConstantDensity::density(Track track_input)
        {
            return constant_density;
        }        
double ConstantDensity::ye(Track track_input)
        {
            return constant_ye;
        }        

/*
----------------------------------------------------------------------
         VariableDensity CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
VariableDensity:VariableDensity(rhomin_input,rhomax_input,ye_input)
        {
            name = "VariableDensity";
            rhomin = rhomin_input;
            rhomax = rhomax_input;
            constant_ye = ye_input;
        }
// track constructor        
VariableDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
VariableDensity::Track::Track()
        {
            //pass
        }
        
double VariableDensity::density(Track track_input)
        {
            // linear density
            m = (rhomax-rhomin)/(track_input.xend-track_input.xini);
            return rhomin+m*(track_input.x-track_input.xini);
        }        
double VariableDensity::ye(Track track_input)
        {
            return constant_ye;
        }
        
/*
----------------------------------------------------------------------
         Star CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
Star:Star(radius_input,iniradius_input,inidensity_input)
        {
            name = "Star";
            constant_ye = 0.86;
            
            radius = radius_input;
            initialradius = iniradius_input;
            initialdensity = inidensity_input;
        }
// track constructor        
Star::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Star::Track::Track()
        {
            //pass
        }
double Star::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Starradius
        double dne;
        double r = radius*x
        if (r <= 1221.50)
            dne = 13.08850-8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150-1.26380*x-3.64260*SQR(x).-5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650-6.47610*x+5.52830*SQR(x)-3.08070*SQR(X)*x.
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970-1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940-8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890-3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910+0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else (r>= radius)
            dne=0.0;
        return dne;
        }
        
double Star::density(Track track_input)
        {
            double rkm = track_input.x/param.km;
            return rdensity(rkm/radius);
        }        
double Star::ye(Track track_input)
        {
            return constant_ye;
        }             
        
/*
----------------------------------------------------------------------
         Earth CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
Earth:Earth()
        {
            name = "Earth";
            radius = 6371.0; // [km]
            constant_ye = 0.494;
        }
// track constructor        
Earth::Track::Track(double xini_input, double xend_input,double baseline_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            baseline = baseline_input;
        }
Earth::Track::Track()
        {
            //pass
        }
double Earth::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        double dne;
        double r = radius*x
        if (r <= 1221.50)
            dne = 13.08850-8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150-1.26380*x-3.64260*SQR(x).-5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650-6.47610*x+5.52830*SQR(x)-3.08070*SQR(X)*x.
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970-1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940-8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890-3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910+0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else (r>= radius)
            dne=0.0;
        return dne;
        }
        
double Earth::density(Track track_input)
        {
            double xkm = track_input.x/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track.baseline/param.km)*xkm);
            return rdensity(r/radius);
        }        
double Earth::ye(Track track_input)
        {
            return constant_ye;
        }             
        
    
        
        
