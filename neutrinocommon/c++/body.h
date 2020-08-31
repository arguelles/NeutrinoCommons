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
        Vacuum(string,double);
        Vacuum();
        
        class Track {
            public :
                double x;
                double xini;
                double xend;                
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
        
        ConstantDensity();
        
        class Track {
            public :
                double x;
                double xini;
                double xend;                
                Track(double,double);
                Track();            
        };
};

class VariableDensity: public Body{
    public:
        double rhomin;
        double rhomax;
        double constant_ye;
        
        VariableDensity();    
        
        class Track {
            public :
                double x;
                double xini;
                double xend;                
                Track(double,double);
                Track();            
        };
};

class Star: public Body{
    public:
        double radius;
        double initialradius;
        double initialdensity;
        double constant_ye;
        
        Star(double,double,double);
        
        double rdensity(double);
        
        class Track {
            public :
                double x;
                double xini;
                double xend;                
                Track(double,double);
                Track();            
        };
};

class Earth: public Body{
    public:
        double radius;
        double constant_ye;        
        Earth();
        
        double rdensity(double);
        
        class Track {
            public :
                double x;
                double xini;
                double xend;
                double baseline;
                Track(double,double,double);
                Track();            
        };
};

// type defining
typedef Body::Track Track;

#endif
