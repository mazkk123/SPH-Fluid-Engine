#ifndef TANK_H_
#define TANK_H_

#include <ngl/VAOFactory.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/Mat4.h>
#include <ngl/Vec3.h>
#include "Globals.h"
#include <vector>

class Tank
{
    public:
        Tank();
        Tank(ngl::Vec3 _p, ngl::Vec3 _d);
        ~Tank();
        //void draw(Camera _c, ngl::Mat4 _mouseTX);
        void buildVAO();
        ngl::Vec3 getPos();
        ngl::Vec3 getDim();
        void setPos();
        void setDim();
        
        //ngl::VAOPrimitives m_t;
        ngl::Vec3 m_TankPos = {0.0f,2.0f,0.0f};
        ngl::Vec3 m_TankDim = {0.5f, 0.5f, 0.5f};
    private:
        Globals globs;
};

#endif