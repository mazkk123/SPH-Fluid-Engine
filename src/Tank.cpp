#include <vector>
#include <ngl/Mat4.h>
#include "Globals.h"
#include "Tank.h"

Tank::Tank()
{

}

Tank::Tank(ngl::Vec3 _p, ngl::Vec3 _d) :
m_TankPos(_p), m_TankDim(_d)
{

}

Tank::~Tank()
{

}

/* void Tank::draw(Camera _c, ngl::Mat4 _mouseTX)
{

} */

void Tank::buildVAO()
{

}

void Tank::setPos()
{
    m_TankPos = globs.tankPos;
}

void Tank::setDim()
{
    m_TankDim = globs.tankDims;
}

ngl::Vec3 Tank::getPos()
{   
    return m_TankPos;
}

ngl::Vec3 Tank::getDim()
{
    return m_TankDim;
}