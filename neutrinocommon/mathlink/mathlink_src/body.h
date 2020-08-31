#ifndef __BODY_H
#define __BODY_H

#include <string>

using namespace std;

class Vacuum {
    public:
        string name;
        double radius;
        
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

// type defining
typedef Vacuum Body;
typedef Vacuum::Track Track;

#endif
