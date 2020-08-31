#ifndef __BODY_H
#define __BODY_H

#include <string>

using namespace std;

class Body{
    public:
        string name;
        
        Body();
        
        class Track{
            public:
                double x;
                double xini;
                double xend;
                Track(double,double);
                Track();
        };
        
        double density(Track);
        double ye(Track);
};

class Vacuum: public Body {
    public:        
        Vacuum();
        
        class Track: public Body::Track {
            public :      
                Track(double,double);
                Track();
        };
        
        double density(Track);
        double ye(Track);
};

class ConstantDensity: public Body{
    public:
        double constant_density;
        double constant_ye;
        
        ConstantDensity(double,double);
        
        class Track: public Body::Track{
            public :              
                Track(double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);
};

class VariableDensity: public Body{
    public:
        double rhomin;
        double rhomax;
        double constant_ye;
        
        VariableDensity(double,double,double);    
        
        class Track: public Body::Track{
            public :          
                Track(double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);        
};

class Star: public Body{
    public:
        double radius;
        double initialradius;
        double initialdensity;
        double constant_ye;
        
        Star(double,double,double);
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :              
                Track(double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);  
};

class Earth: public Body{
    public:
        double radius;
        double constant_ye;        
        Earth();
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :
                double baseline;
                Track(double,double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);          
};

/*
int GetTransversedWidthAndDistance(double[],double,double);

class EarthWithCavity: public Body{
    public:
        double radius;
        double constant_ye;
        // cavity parameters description
        // cavity matter properties
        double rho_cavity;
        double ye_cavity;
        // cavity geometry
        double wx;
        double wy;
        double th_rot;
        double deep_cavity;
        double th_cavity;        
        
        EarthWithCavity(double, double, double, double, double, double, double);
        
        double rdensity(double);
        
        class Track {
            public :
                double phi;
                double th;
                double baseline_angle;
                
                double x;
                double xini;
                double xend;
                double baseline;
                
                double xin;
                double xout;
                Track(EarthWithCavity,double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);
        bool isInCavity(Track);
        double GetBaseline(Track);
        int GetTransversedWidthAndDistance(double[],Track);
};
*/

class AGN: public Body{
    public:
        double radius;
        double initialdensity;
        double constant_ye;
        
        AGN(double,double,double);
        
        double rdensity(double);
        
        class Track: public Body::Track{
            public :            
                Track(double,double);
                Track();            
        };
        
        double density(Track);
        double ye(Track);  
};

// type defining
typedef Body::Track Track;

#endif
