#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string_view>
#include <vector>
#include <ngl/Vec3.h>

class Globals
{
    public:
        Globals();
        ~Globals();

        float dt = 0.02f;
        ngl::Vec3 tankDims = {3.0f, 3.0f, 3.0f};
        ngl::Vec3 tankPos = {0.0f , 2.0f, 0.0f};

        unsigned int total_particles = 2300;
        float mass = 0.02f;

        unsigned int p1 = 73856093;
        unsigned int p2 = 19349663;
        unsigned int p3 = 83492791;

        float _rd = 998.2f;
        float _kconst = 5;
        float _mvisc = 3.5f;
        float _sigmasurf = 0.0728f;
        float _lsurf = 5.0f;
        float _buoyancyConst = 0.0f;
        float _h = 0.11f;
        int particleSize = 2;
        float particleSep = 0.012f;
        float RLOS = 0.03f;
        ngl::Vec3 _dims = {15.0f, 15.0f ,15.0f};
};

#endif




