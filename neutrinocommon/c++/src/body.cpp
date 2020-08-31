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
ConstantDensity::ConstantDensity(double density_input,double ye_input)
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
VariableDensity::VariableDensity(double rhomin_input,double rhomax_input,double ye_input)
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
            double m;
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
Star::Star(double radius_input,double iniradius_input,double inidensity_input)
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
        // x is adimentional radius : x = 0 : center, x = 1 : Starradius
            double r = radius*x;
            if (r <= radius and r> initialradius){
                return initialdensity*exp( - r/initialradius );
            } else if ( r < initialradius){
                return initialdensity;
            } else
                return 0.0;
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
Earth::Earth()
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
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }
        
double Earth::density(Track track_input)
        {
            double xkm = track_input.x/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track_input.baseline/param.km)*xkm);
            return rdensity(r/radius);
        }        
double Earth::ye(Track track_input)
        {
            return constant_ye;
        }       
        


        
/*
----------------------------------------------------------------------
         EarthWithCavity CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

/*
// auxiliary function

int GetTransversedWidthAndDistance(double wd[],double phi, double th,\
                                   double wx,double wy,double th_rot,double radius,double deep_cavity, double th_cavity)
        {
            //double wd[2];
            wd[0] = 0.0;
            wd[1] = 0.0;
            
            if ( wx != 0.0 && wy != 0.0){
                
                double rth = th_rot;

                double xs = -radius*sin(phi);
                double ys = radius*cos(phi);
                
                double R_cavity = (radius - deep_cavity);
                
                double xc = -R_cavity*tan(th_cavity);
                double yc = R_cavity;
                
                double wxt = wx/2.0;
                double wyt = wy/2.0;
                
                double m = tan(th);
                
                double cxx = SQR(cos(rth))/SQR(wxt) + SQR(sin(rth))/SQR(wyt);
                double cyy = SQR(sin(rth))/SQR(wxt) + SQR(cos(rth))/SQR(wyt);
                double cxy = sin(2.0*rth)*(1.0/(wxt*wyt) - 1.0/(wxt*wxt));
                
                double cx = -(2.0*cxx*xc+cxy*yc);
                double cy = -(2.0*cyy*yc+cxy*xc);
                double cc = cxx*SQR(xc)+cyy*SQR(yc)+xc*yc*cxy-1.0;
                
                double a2 = m*m*cyy+m*cxy+cxx;
                double b2 = 2*m*(ys-xs*m)*cyy+(ys-xs*m)*cxy+m*cy+cx;
                double c2 = SQR(ys-xs*m)*cyy+cy*(ys-xs*m)+cc;
                double delta2 = b2*b2-4*a2*c2;
                
                if (delta2 <= 0.0){
                    return 1;
                } else {
                double c1x = -b2/(2.0*a2)+sqrt(delta2)/(2.0*a2);
                double c1y = ys + m*(c1x-xs);
                double c2x = -b2/(2.0*a2)-sqrt(delta2)/(2.0*a2);
                double c2y = ys + m*(c2x-xs);

                double w2 = sqrt(SQR(c1x-c2x)+SQR(c1y-c2y));

                double d2 = sqrt(SQR(c2x-xs)+SQR(c2y-ys));
                
                wd[0] = w2/param.km;
                wd[1] = d2/param.km;
                return 0;
                }
                
            } else {
                return 1;
            }
        }                

// constructor
EarthWithCavity::EarthWithCavity(double rho_cavity_input,double ye_cavity_input,double wx_input,double wy_input,double th_rot_input,double deep_cavity_input,double th_cavity_input)
        {
            name = "EarthWithCavity";
            radius = 6371.0*param.km; // [km]
            constant_ye = 0.494;
            rho_cavity = rho_cavity_input;
            ye_cavity = ye_cavity_input;
            wx = wx_input;
            wy = wy_input;
            th_rot = th_rot_input;
            deep_cavity = deep_cavity_input;
            th_cavity = th_cavity_input;
        }
// track constructor        
EarthWithCavity::Track::Track(EarthWithCavity earth,double phi_input, double th_input)
        {
            double radius = 6371.0*param.km; // [km]
            phi = phi_input;
            th = th_input;
            baseline_angle = 2.0*(phi-th);
            double L = sqrt(2.0*SQR(radius)*(1.0-cos(baseline_angle)));
            
            x = 0.0;
            xini = 0.0;
            xend = L;
            
            double wd[2];
            GetTransversedWidthAndDistance(wd,phi,th,\
                                           earth.wx,earth.wy,earth.th_rot,earth.radius,earth.deep_cavity,earth.th_cavity);
            xini = wd[0]*param.km;
            xout = (wd[0]+wd[1])*param.km;
        }
EarthWithCavity::Track::Track()
        {
            //pass
        }
double EarthWithCavity::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : EarthWithCavityradius
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }
        
double EarthWithCavity::density(Track track_input)
        {
            if (track_input.x >= track_input.xin && track_input.x <= track_input.xout){
                return rho_cavity;    
            } else {
                double xkm = track_input.x/param.km;
                double r = sqrt(SQR(radius)+SQR(xkm)-(track_input.baseline/param.km)*xkm);
                return rdensity(r/radius);
            }

        }        
double EarthWithCavity::ye(Track track_input)
        {
            if (track_input.x >= track_input.xin && track_input.x <= track_input.xout){
                return ye_cavity;    
            } else {
                return 0.494;
            }
        }
        
bool EarthWithCavity::isInCavity(Track track_input)
        {
            double phi = track_input.phi;
            double th = track_input.th;
            
            double x_source = -radius*sin(phi);
            double y_source = radius*cos(phi);
            
            double xcurrent = track_input.x;
            double x_neutrino = x_source + xcurrent*cos(th);
            double y_neutrino = y_source + xcurrent*sin(th);
            
            double R_cavity = (radius - deep_cavity);
                
            double x_cavity = -R_cavity*tan(th_cavity);
            double y_cavity = R_cavity;
            
            double x_rot = x_neutrino*cos(th_rot) - y_neutrino*sin(th_rot);
            double y_rot = x_neutrino*sin(th_rot) + y_neutrino*cos(th_rot);
            
            double a_cav = wx/2.0;
            double b_cav = wy/2.0;
            
            return SQR((x_rot-x_cavity)/a_cav)+SQR((y_rot-y_cavity)/b_cav) <= 1.0;
        }
        
double EarthWithCavity::GetBaseline(Track track_input)
        {
            double baseline_angle = 2.0*(track_input.phi-track_input.th);
            return sqrt(2.0*SQR(radius)*(1.0 - cos(baseline_angle)));
        }
        
int EarthWithCavity::GetTransversedWidthAndDistance(double wd[],Track track_input)
        {
            //double wd[2];
            wd[0] = 0.0;
            wd[1] = 0.0;
            
            if ( wx != 0.0 && wy != 0.0){
                
                double phi = track_input.phi;
                double th = track_input.th;
                double rth = th_rot;

                double xs = -radius*sin(phi);
                double ys = radius*cos(phi);
                
                double R_cavity = (radius - deep_cavity);
                
                double xc = -R_cavity*tan(th_cavity);
                double yc = R_cavity;
                
                double wxt = wx/2.0;
                double wyt = wy/2.0;
                
                double m = tan(th);
                
                double cxx = SQR(cos(rth))/SQR(wxt) + SQR(sin(rth))/SQR(wyt);
                double cyy = SQR(sin(rth))/SQR(wxt) + SQR(cos(rth))/SQR(wyt);
                double cxy = sin(2.0*rth)*(1.0/(wxt*wyt) - 1.0/(wxt*wxt));
                
                double cx = -(2.0*cxx*xc+cxy*yc);
                double cy = -(2.0*cyy*yc+cxy*xc);
                double cc = cxx*SQR(xc)+cyy*SQR(yc)+xc*yc*cxy-1.0;
                
                double a2 = m*m*cyy+m*cxy+cxx;
                double b2 = 2*m*(ys-xs*m)*cyy+(ys-xs*m)*cxy+m*cy+cx;
                double c2 = SQR(ys-xs*m)*cyy+cy*(ys-xs*m)+cc;
                double delta2 = b2*b2-4*a2*c2;
                
                if (delta2 <= 0.0){
                    return 1;
                } else {
                double c1x = -b2/(2.0*a2)+sqrt(delta2)/(2.0*a2);
                double c1y = ys + m*(c1x-xs);
                double c2x = -b2/(2.0*a2)-sqrt(delta2)/(2.0*a2);
                double c2y = ys + m*(c2x-xs);

                double w2 = sqrt(SQR(c1x-c2x)+SQR(c1y-c2y));

                double d2 = sqrt(SQR(c2x-xs)+SQR(c2y-ys));
                
                wd[0] = w2/param.km;
                wd[1] = d2/param.km;
                return 0;
                }
                
            } else {
                return 1;
            }
        }
        
*/        
        
/*
----------------------------------------------------------------------
         AGN CLASS DEFINITIONS
----------------------------------------------------------------------         
*/

// constructor
AGN::AGN(double radius_input,double inidensity_input, double ye_input)
        {
            name = "AGN";
                        
            radius = radius_input;
            initialdensity = inidensity_input;
            constant_ye = ye_input;            
        }
// track constructor        
AGN::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
AGN::Track::Track()
        {
            //pass
        }
double AGN::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : AGNradius
            return initialdensity*exp( -x * 10);
        }
        
double AGN::density(Track track_input)
        {
            double rparsec = track_input.x/param.parsec;
            return rdensity(rparsec/radius);
        }        
double AGN::ye(Track track_input)
        {
            return constant_ye;
        }        