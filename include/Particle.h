#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <ngl/Mat4.h>
#include <ngl/Vec3.h>
#include <array>
#include "Globals.h"
#include "Tank.h"
#include <iostream>

class Particle
{
    public:
        Particle() = default;
        ~Particle();

        void draw(ngl::Mat4 _mouseTX);
        void setId(const unsigned int &_i);
        unsigned int &getId();
        void setHashKey(const unsigned int &_hashKey);
        void setPrevHashKey(const unsigned int &_prevHashKey);
        unsigned int getPrevHashKey();
        unsigned int &getHashKey();
        std::vector<unsigned int>& getNeighbours();
        void setNeighbours(const std::vector<unsigned int> &_neighbours);

        void setLife(const unsigned int &_life);
        unsigned int getLife();
        void setMaxLife(const unsigned int &_maxLife);
        unsigned int getMaxLife();
        void randomVectorOnSphere(float _radius );

        void setPosition(const ngl::Vec3 &_pos);
        ngl::Vec3 &getPosition();
        void setMass();
        float &getMass();
        void setVelocity(const ngl::Vec3 &_vel);
        void setVelocity(const float &_vel);
        ngl::Vec3 &getVelocity();
        void setXSPHVelocity();
        void setAcceleration(const ngl::Vec3 &_acc);
        ngl::Vec3 &getAcceleration();
        void setDensity(const float &_md);
        float &getDensity();
        float &getPressure();
        void setPressure(const float &_pressure);
        ngl::Vec3 &getPressureForce();
        void setPressureForce(const ngl::Vec3 &_fp);
        void setViscosity(const ngl::Vec3 &_fv);
        ngl::Vec3 &getViscosity();
        void setGravityForce(const ngl::Vec3 &_fg);
        ngl::Vec3 &getGravityForce();
        void setBuoyancyForce(const ngl::Vec3 &_fb);
        ngl::Vec3 &getSumForces();
        void setSumForces(const ngl::Vec3 &_sumForces);
        ngl::Vec3 &getBuoyancyForce();
        void setNormalField(const ngl::Vec3 &_normalField);
        ngl::Vec3 &getNormalField();
        void setLaplNormalField(const float &_normalLaplField);
        float &getLaplNormalField();
        void setSurfaceCurvature(const float &_sufaceCurvature);
        float &getSurfaceCurvature();
        void setSurfaceTension(const ngl::Vec3 &_fst);
        ngl::Vec3 &getSurfaceTension();
        void setXSPHVel(const ngl::Vec3 &_XSPHVel);
        ngl::Vec3 &getXSPHVel();

        float mass;
        float m_md;
        float m_pressure;

        ngl::Vec3 m_pos;
        ngl::Vec3 m_acc;
        ngl::Vec3 m_vel;
        const ngl::Vec3 m_colour = {0.0f, 0.0f, 0.0f};

        ngl::Vec3 m_fp;
        ngl::Vec3 m_fv;
        ngl::Vec3 m_fg = {0.0f, 1.2f, 0.0f};
        ngl::Vec3 m_fb;
        ngl::Vec3 m_fst;
        ngl::Vec3 m_sumForces;
        ngl::Vec3 m_XSPHVel;
        ngl::Vec3 colour;
        
	    unsigned int m_life = 0;
	    unsigned int m_maxLife = 0;
        unsigned int m_id;
        unsigned int m_hashKey;
        unsigned int m_prevHashKey;

        float m_normalLaplField;
        float m_surfaceCurvature;
        ngl::Vec3 m_normalField;

        std::vector<unsigned int>m_neighbours={};

        Globals globs;

};

#endif
