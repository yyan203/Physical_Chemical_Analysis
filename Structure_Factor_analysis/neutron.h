//
// Created by yxy277 on 3/21/18.
//

#include <map>

class Neutron {
private:
    std::map< int,double > b; // Neutron scattering lengths, from www.ncnr.nist.gov/resources/n-lengths

public:
    Neutron (){
        /* debug */
        //b[1] = 1.0 ; // C
        //b[2] = 1.0 ; // D,i.e.2H yes, it is negative
        //b[3] = 1.0 ; // O
        //b[4] = 1.0 ; // N
        //b[5] = 1.0 ; // Zn

        ////////////////////////////////  for Si-O system
        b[1] = 5.803 ; // O
        //b[2] = -3.739; // H, yes, it is negative
        b[2] = 4.1491; // Si
        // b[3] = 5.803 ; // O
        // b[4] = 9.36  ; // N
        // b[5] = 5.680 ; // Zn
    }
    double getb(int i) {
        if (i<1 || i > 5) {printf("index error!"); exit(0);}
        return b[i];
    }
    };

