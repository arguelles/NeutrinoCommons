#include "body.h"

// vacuum constructor
Vacuum::Vacuum(string name_input = "Vacuum", double radius_input = 1.0)
        {
            name = name_input;
            radius = radius_input;
        }
Vacuum::Vacuum()
        {
            //pass
            name = "Vacuum";
            radius = 1.0;
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
