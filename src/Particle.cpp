#include <vector>
#include "System.h"
#include "Globals.h"
#include "Particle.h"
#include <math.h>


Particle::~Particle()
{
  // Particle destructor
}

//CALCULATES THE PHYSICAL FORCES TO APPLY TO PARTICLES

void Particle::setId(const unsigned int &_i)
{
    m_id = _i;
}

unsigned int &Particle::getId()
{
    return m_id;
}

unsigned int &Particle::getHashKey()
{
    return m_hashKey;
}

unsigned int Particle::getPrevHashKey()
{
    return m_prevHashKey;
}

void Particle::setPrevHashKey(const unsigned int &_prevHashKey)
{
    m_prevHashKey = _prevHashKey;
}

void Particle::setHashKey(const unsigned int &_hashKey)
{
    m_hashKey = _hashKey;
}

void Particle::setLife(const unsigned int &_life)
{
    m_life = _life;
}

unsigned int Particle::getLife()
{
    return m_life;
}

void Particle::setMaxLife(const unsigned int &_maxLife)
{
    m_maxLife = _maxLife;
}

unsigned int Particle::getMaxLife()
{
    return m_maxLife;
}

std::vector<unsigned int>& Particle::getNeighbours()
{
    return m_neighbours;
}

void Particle::setNeighbours(const std::vector<unsigned int> &_neighbours)
{
    for(auto &i: _neighbours)
    {
        m_neighbours.push_back(i);
    }
}

void Particle::setPosition(const ngl::Vec3 &_pos)
{
    m_pos = _pos;
}

ngl::Vec3& Particle::getPosition()
{
    return m_pos;
}

void Particle::setMass()
{
    mass = globs.mass;
}

float& Particle::getMass()
{
    return mass;
}
 
void Particle::setVelocity(const ngl::Vec3 &_vel)
{
    m_vel = _vel;
}

void Particle::setVelocity(const float &_vel)
{
    m_vel = {_vel, _vel, _vel};
}

ngl::Vec3& Particle::getVelocity() 
{
    return m_vel;
}

void Particle::setAcceleration(const ngl::Vec3 &_acc)
{
    m_acc = _acc;
}

ngl::Vec3& Particle::getAcceleration() 
{
    return m_acc;
}

void Particle::setDensity(const float &_md)
{
    m_md = _md;
}   

float& Particle::getDensity() 
{
    return m_md;
}

float &Particle::getPressure()
{
    return m_pressure;
}

void Particle::setPressure(const float &_pressure)
{
    m_pressure = _pressure;
}

void Particle::setPressureForce(const ngl::Vec3 &_fp)
{
    m_fp = _fp;
}

ngl::Vec3 &Particle::getPressureForce()
{
    return m_fp;
}

void Particle::setViscosity(const ngl::Vec3 &_fv)
{
    m_fv = _fv;
}

ngl::Vec3 &Particle::getViscosity()
{
    return m_fv;
}

void Particle::setGravityForce(const ngl::Vec3 &_fg)
{
    m_fg = _fg;
}

ngl::Vec3 &Particle::getGravityForce()
{
    return m_fg;
}

void Particle::setBuoyancyForce(const ngl::Vec3 &_fb)
{
    m_fb = _fb;
}

ngl::Vec3 &Particle::getBuoyancyForce()
{
    return m_fb;
}

void Particle::setNormalField(const ngl::Vec3 &_normalField)
{
    m_normalField = _normalField;
}

ngl::Vec3 &Particle::getNormalField()
{
    return m_normalField;
}

float &Particle::getLaplNormalField()
{
    return m_normalLaplField;
}

void Particle::setLaplNormalField(const float &_normalLaplField)
{
    m_normalLaplField = _normalLaplField;
}

void Particle::setSurfaceCurvature(const float &_sufaceCurvature)
{
    m_surfaceCurvature = _sufaceCurvature;
}

float &Particle::getSurfaceCurvature()
{
    return m_surfaceCurvature;
}

void Particle::setSurfaceTension(const ngl::Vec3 &_fst)
{
    m_fst = _fst;
}

ngl::Vec3 &Particle::getSurfaceTension()
{
    return m_fst;
}

void Particle::setSumForces(const ngl::Vec3 &_sumForces)
{
    m_sumForces = _sumForces;
}

ngl::Vec3 &Particle::getSumForces()
{
    return m_sumForces;
}

void Particle::setXSPHVel(const ngl::Vec3 &_XSPHVel)
{
    m_XSPHVel = _XSPHVel;
}

ngl::Vec3 &Particle::getXSPHVel()
{
    return m_XSPHVel;
}
