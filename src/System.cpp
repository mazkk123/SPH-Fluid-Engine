#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include "Globals.h"
#include "System.h"
#include "Particle.h"
#include <ngl/Random.h>
#include <ngl/NGLStream.h>
#include <ngl/SimpleVAO.h>
#include <ngl/VAOFactory.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <fmt/format.h>

//Update Particle Operations

System::System( ngl::Vec3 _pos, unsigned int frame)
{
    m_pos = _pos;
	m_numParticles = globs.total_particles;
    m_frame = frame;

    // zero the position of every particle in the system
    // set the number of particles in the system to the global
    // particle amount
    // set the number of frames to render as the input frame number

    particleReset(m_numParticles);
    // calls a separate function to reset all particles into a tank
    // and space them into a grid of uniformly spaced particles with 
    // dimensions specified by the user

    m_vao = ngl::VAOFactory::createVAO(ngl::simpleVAO, GL_POINTS);
    // creates a vertex array object to be used in a subsequent operation
    // for drawing particles onto the screen
}

void System::particleReset(unsigned int &_numParticles)
{
    // on initialization and user restart system calls, resets
    // particle positions to be within a uniformly spaced grid of 
    // particles, sets their id and adds particle ids to their 
    // appropriate cells in the hash table

    m_particles.resize(_numParticles);
    resetParticleToTank();

    unsigned int iterator = 0;
    for(auto &p: m_particles)
    {
        p.setId(iterator);
        ids.push_back(p.getId());
        iterator += 1;
    }
    addToHashTable();
}

void System::findNeighbours(const unsigned int &_id)
{
    // finds particle neighbours from a given particle id and 
    // sets the particle neighbour list with particle neighbours sourced 
    // from the hash map of particle ids within the same cell

    std::vector<unsigned int>_neighbours={};
    _neighbours.clear();
    for(auto &data: getHashTable()[getParticle(_id).getHashKey()])
    {
        for(auto &id: ids)
        {
            if(data == id)
            {
                _neighbours.push_back(data);
                break;
            }
        }
    }

    if (!_neighbours.empty()) // Check if the vector is not empty
    {
        _neighbours.erase(_neighbours.begin());
    }

    getParticle(_id).setNeighbours(_neighbours);
}

System::~System()
{
  // class destructor
}

void System::update()
{
    // called on update every frame
    // iterates through every particle in the std::vector m_particles
    // updates their hash table key, sets their mass and updates their position
    for (auto &p : m_particles)
    {
        findNeighbours(p.getId());
        p.setMass();

        updatePosition(p.getId());
        updateHashTable(p);
    }
}

void System::resetParticleToTank()
{

    unsigned int rows = globs._dims.m_x;
    unsigned int columns = globs._dims.m_y;
    unsigned int depth = globs._dims.m_z;
    float _particleSep = globs.particleSep;
    unsigned int particleIndex = 0;

    for (float row = 0.0f; row<rows; ++row)
    {
        for (float column = 0.0f; column<columns; ++column)
        {
            for (float d = 0.0f; d<depth; ++d)
            {
                if (particleIndex < m_particles.size())
                {
                    ngl::Vec3 _tempPos = {row*_particleSep, 
                                          column*_particleSep, 
                                          d*_particleSep};
                    m_particles[particleIndex].setPosition(_tempPos);
                    particleIndex++;
                }
                else
                {
                    break;
                }
            }
        }
    }
    
}


void System::resetParticleToSphere()
{

}

void System::draw()
{
    // create a new flat array for the particles
    if(m_particles.size()==0)
        return;

    m_vao->bind();
	m_vao->setData(ngl::SimpleVAO::VertexData(m_particles.size()*sizeof(Particle), m_particles[0].m_pos.m_x));
	m_vao->setVertexAttributePointer(0, 3, GL_FLOAT,sizeof(Particle), 0);
	m_vao->setNumIndices(m_particles.size());
	m_vao->draw();
	m_vao->unbind();
}

ngl::Vec3 System::randomVectorOnSphere(float _radius)
{
	
	float phi = ngl::Random::randomPositiveNumber(static_cast<float>(M_PI * 2.0f));
	float costheta = ngl::Random::randomNumber();
	float u = ngl::Random::randomPositiveNumber();
	float theta = acos(costheta);
	float r = _radius * std::cbrt(u);
	return ngl::Vec3(r * sin(theta) * cos(phi),
					 r * sin(theta) * sin(phi),
					 r * cos(theta) );
}

//MATHS CALCULATIONS FOR SPATIAL HASHING SEARCH
bool System::isPrime(unsigned int &_numToCheck)
{
    if (_numToCheck<=1)
    {
        return false;
    }
    for (int i = 2; i <= _numToCheck/2; ++i)
    {
        if (_numToCheck % i == 0)
        {
            return false;
        }
    }
    return true;
}

unsigned int System::primeFunction(unsigned int _n)
{
    // finds the next prime number after the input number _n
    unsigned int nextNum = _n + 1;
    while(!isPrime(nextNum))
    {
        nextNum++;
    }    
    return nextNum;
}

ngl::Vec3 System::x_hat(const ngl::Vec3 &m_pos)
{
    // applies the hash function on the input ngl::Vec3 particle
    // position

    ngl::Vec3 local_x_hat = {std::floor(m_pos[0]/globs._h), 
                             std::floor(m_pos[1]/globs._h),
                             std::floor(m_pos[2]/globs._h)};
    return local_x_hat;
}

//HASH TABLE METHODS
void System::addToHashTable()
{   
    // initially adds all particles to their respective cell in the hash map
    // by calculating their hash key, and finding whether a cell with that hash key
    // already exists in the hash map. If it does, the particle id is added to that cell
    // otherwise, it is added to a new cell if one doesn't already exist with a matching
    // hash key

    for(auto &p: m_particles)
    {
        ngl::Vec3 local_x_hat = x_hat(p.getPosition());
        unsigned int x_hat_0 = (static_cast<int>(local_x_hat[0])*globs.p1)^
                               (static_cast<int>(local_x_hat[1])*globs.p2)^
                               (static_cast<int>(local_x_hat[2])*globs.p3);
        unsigned int _kHash = x_hat_0*primeFunction(2*m_numParticles);
        p.setHashKey(_kHash);

        for(auto &cell: getHashTable())
        {
            if (cell.first == _kHash)
            {
                m_hashMap[cell.first].push_back(p.getId());
            }
        }
        if(m_hashMap[_kHash].begin()==m_hashMap[_kHash].end())
        {
            m_hashMap[_kHash]={p.getId()};
        }       
    }    
}

void System::updateHashTable(Particle &_p)
{
    // called on particle update to replace the hash map cells 
    // based on updated particle positions and hash keys for every frame
    // of the simulation

    ngl::Vec3 local_x_hat = x_hat(_p.getPosition());
    unsigned int x_hat_0 = (static_cast<int>(local_x_hat[0])*globs.p1)^
                           (static_cast<int>(local_x_hat[1])*globs.p2)^
                           (static_cast<int>(local_x_hat[2])*globs.p3);

    unsigned int _kHash = x_hat_0*primeFunction(2*m_numParticles);
    _p.setHashKey(_kHash);

    bool foundMatchingID = false;
    for(auto &cell: getHashTable())
    {
        if (cell.first == _kHash)
        {
            for(auto &cell_val : cell.second)
            {
                if(_p.getId() == cell_val)
                {
                    foundMatchingID = true;
                    break;
                }
            }

            if (!foundMatchingID)
            {
                m_hashMap[cell.first].push_back(_p.getId());
            }
        }
    }
}

//GETTER AND SETTER METHODS 
Particle &System::getParticle(const unsigned int &_i)
{
    // retrieve a particle from m_particles from its id
    return m_particles[_i];
}

std::unordered_map<unsigned int, std::vector<unsigned int>> System::getHashTable()
{
    // retrieve a copy of the hash map
    return m_hashMap;
}

//FINDS OUT THE WEIGHTING FUNCTIONS IN SPH APPROXIMATION.

float System::kernel(const ngl::Vec3 &_r, const int &_type)
{
    float _kernel = 0.0f;
    float _length = _r.length();

    if(_length >= 0.0f && _length <= globs._h)
    {
        if (_length == 0.0f || std::isnan(_length)) {
            return 0.0f;
        }

        switch(_type)
        {
            case 0:
            {
                float _poly6Cutoff = std::pow(globs._h, 2);
                float _poly6Cutoff2 = std::pow(_length, 2);
                float _poly6Dist = std::pow((_poly6Cutoff - _poly6Cutoff2), 3);
                _kernel = (315.0f / (64.0f * M_PI * std::pow(globs._h, 9))) * _poly6Dist;
                _kernel/= (315.0f / (64.0f * M_PI * std::pow(globs._h, 9)));
                break;
            }
            case 1:
            {
                float _spikyCutoff = globs._h;
                float _spikyDist = _spikyCutoff - _length;
                _kernel = (15.0f / (M_PI * std::pow(globs._h, 6))) * std::pow(_spikyDist, 3);
                _kernel /= (15.0f / (M_PI * std::pow(globs._h, 6)));
                break;
            }
            case 2:
            {
                float _viscosityCutoff = globs._h;
                float _viscosityCutoff2 = std::pow(globs._h, 2);
                float _viscosityDist = (_viscosityCutoff - _length) / _viscosityCutoff;
                _kernel = (15.0f / (2.0f * M_PI * std::pow(globs._h, 3))) * (-1 * std::pow(_viscosityDist, 3) + std::pow(_viscosityDist, 2) + (globs._h / (_length * 2.0f)) - (globs._h / (3.0f * _viscosityCutoff2)));
                _kernel /= (15.0f / (2.0f * M_PI * std::pow(globs._h, 3))) ;
                break;
            }
            default:
                break;
        }
    }

    if (std::isnan(_kernel)) {
        return 0.0f;
    }

    // Normalize the kernel value to be between 0 and 1
    _kernel = std::min(_kernel, 1.0f);
    _kernel = std::max(_kernel, 0.0f);

    return _kernel;
}

ngl::Vec3 System::kernelGradient(const ngl::Vec3 &_r, const int &_type)
{
    float _dist = _r.length();

    float _kMultPoly6=-945.0f/(32.0f*M_PI*std::pow(globs._h, 9));
    float _kDistBetweenMultPoly6 = std::pow((std::pow(globs._h,2) - std::pow(_dist,2)),2);

    float _kMultSpiky=-45.0f/(M_PI*std::pow(globs._h, 6));
    float _kDistBetweenMultSpiky = std::pow((globs._h - _dist),2);

    float _kMultViscosity=15.0f/(2.0f*M_PI*std::pow(globs._h, 3));
    float _kDistBetweenMultViscosity = ((-3*_dist/(2*std::pow(globs._h,3))) +
                                            2/std::pow(globs._h, 2)
                                            - globs._h/(2*std::pow(_dist, 3)));

    ngl::Vec3 _kGrad(0, 0, 0);
    switch(_type)
    {
        case 0:
            _kGrad.set(_kMultPoly6*_kDistBetweenMultPoly6*_r[0],
                        _kMultPoly6*_kDistBetweenMultPoly6*_r[1],
                        _kMultPoly6*_kDistBetweenMultPoly6*_r[2]);
            break;
        case 1:
            _kGrad.set(_kMultSpiky*(1/_dist)*_kDistBetweenMultSpiky*_r[0],
                        _kMultSpiky*(1/_dist)*_kDistBetweenMultSpiky*_r[1],
                        _kMultSpiky*(1/_dist)*_kDistBetweenMultSpiky*_r[2]);
            break;
        case 2:
            _kGrad.set(_kMultViscosity*_kDistBetweenMultViscosity*_r[0],
                        _kMultViscosity*_kDistBetweenMultViscosity*_r[1],
                        _kMultViscosity*_kDistBetweenMultViscosity*_r[2]);
            break;
        default:
            break;
    }

    // Normalize the kernel gradient
    //_kGrad.normalize();

    return _kGrad;
}


float System::kernelLaplacian(const ngl::Vec3 &_r, const int &_type)
{
    float _kernelLapl;

    float _kernelLaplPoly6M = -945.0f/(32*M_PI*std::pow(globs._h,9));
    float _kernelLaplPoly6DistBetw = (std::pow(globs._h,2) - std::pow(_r.length(),2))*
                                     (3*std::pow(globs._h,2) - 7*std::pow(_r.length(),2));
    
    float _kernelLaplSpiky = -90.0f/(M_PI*std::pow(globs._h,6));
    float _kernelLaplSpikyDistBetw = (globs._h - _r.length()) * (globs._h - 2*_r.length());
    
    float _kernelLaplVisc = 45.0f/(M_PI * std::pow(globs._h,6));
    float _kernelLaplViscDistBetw = (globs._h - _r.length());

    switch(_type)
    {
        case 0:
            _kernelLapl = _kernelLaplPoly6M*_kernelLaplPoly6DistBetw;
            break;
        case 1: 
            _kernelLapl = _kernelLaplSpiky*_kernelLaplSpikyDistBetw;
            break;
        case 2:
            _kernelLapl = _kernelLaplVisc*_kernelLaplViscDistBetw;
            break;
    }
    
    return _kernelLapl;
}

//PARTICLE UPDATES

void System::updatePosition(const unsigned int &_id)
{
    updateVelocity(_id);
    TankCollision(_id);
    updateXSPHVelocity(_id);
    ngl::Vec3 _position = getParticle(_id).getPosition() + globs.dt*getParticle(_id).getVelocity();
    getParticle(_id).setPosition(_position);
}

void System::updateVelocity(const unsigned int &_id)
{
    updateAcceleration(_id);
    ngl::Vec3 _vel = getParticle(_id).getVelocity() + globs.dt* getParticle(_id).getAcceleration();
    getParticle(_id).setVelocity(_vel);
}

void System::updateXSPHVelocity(const unsigned int &_id)
{
    float _XSPHVel = 0;
    for(const auto &i: getParticle(_id).getNeighbours())
    {
        float _neighbourMassDensity = getParticle(i).getDensity();
        float _neighbourMass = getParticle(i).getMass();
        ngl::Vec3 _neighbourPos = getParticle(i).getPosition();
        ngl::Vec3 _posDifference = getParticle(_id).getPosition() - _neighbourPos;
        _XSPHVel += 2.0f*_neighbourMass / (getParticle(_id).getDensity() + _neighbourMassDensity) * kernel(_posDifference, 0);
    }
    ngl::Vec3 _XSPHFinal = {getParticle(_id).getVelocity()[0] + globs._sigmasurf * _XSPHVel,
                            getParticle(_id).getVelocity()[1] + globs._sigmasurf * _XSPHVel,
                            getParticle(_id).getVelocity()[2] + globs._sigmasurf * _XSPHVel};   
    getParticle(_id).setXSPHVel(_XSPHFinal);
}

void System::updateAcceleration(const unsigned int &_id)
{
    updateSumForces(_id);
    ngl::Vec3 _acc = getParticle(_id).getSumForces()/getParticle(_id).getDensity();
    getParticle(_id).setAcceleration(_acc);
}

void System::updateSumForces(const unsigned int &_id)
{
    updateDensity(_id);
    updateBuoyancyForce(_id);
    updatePressureForce(_id);
    updateSurfaceTension(_id);
    updateViscosity(_id);
    updateGravityForce(_id);
    
    ngl::Vec3 _surfTension = getParticle(_id).getSurfaceTension();
    ngl::Vec3 _viscosity = getParticle(_id).getViscosity();
    ngl::Vec3 _buoyancy = getParticle(_id).getBuoyancyForce();
    ngl::Vec3 _pressureForce = getParticle(_id).getPressureForce();
    ngl::Vec3 _gravityForce = getParticle(_id).getGravityForce();
    ngl::Vec3 _sumForces = _gravityForce + _buoyancy + _surfTension + _viscosity + _pressureForce;
    getParticle(_id).setSumForces(_sumForces);
}


void System::updateDensity(const unsigned int &_id)
{
    float _densityFinal = 0;
    float _neighbourMass = 0;
    ngl::Vec3 _posDifference = {0,0,0};
    ngl::Vec3 _neighbourPos = {0,0,0};
    float k = 0;

    for(const auto &i: getParticle(_id).getNeighbours())
    { 
        _neighbourMass = getParticle(i).getMass();
        _neighbourPos = getParticle(i).getPosition();
        _posDifference = getParticle(_id).getPosition() - _neighbourPos;
        k = kernel(_posDifference, 0);
        _densityFinal += _neighbourMass*k;
    }
    _densityFinal += globs._rd;
    getParticle(_id).setDensity(_densityFinal);
}   


void System::updatePressure(const unsigned int &_id)
{
    float _pressure = (globs._kconst*getParticle(_id).getDensity()) - globs._kconst*globs._rd;
    getParticle(_id).setPressure(_pressure);
}

void System::updatePressureForce(const unsigned int &_id)
{
    updatePressure(_id);
    ngl::Vec3 _pressureForceFinal = {0,0,0};
    float _neighbourPressure = 0;
    float _neighbourMassDensity = 0;
    float _neighbourMass = 0;
    ngl::Vec3 _neighbourPos = {0,0,0};
    ngl::Vec3 _posDifference = {0,0,0};
    ngl::Vec3 _gradPosDifference = {0,0,0};

    for(const auto &i: getParticle(_id).getNeighbours())
    {
        _neighbourPressure = getParticle(i).getPressure();
        _neighbourMassDensity = getParticle(i).getDensity();
        _neighbourMass = getParticle(i).getMass();
        _neighbourPos = getParticle(i).getPosition();
        _posDifference = getParticle(_id).getPosition() - _neighbourPos;
        _gradPosDifference = kernelGradient( _posDifference, 1);
        
        ngl::Vec3 _pressureForceTemp;

        if(std::isnan(_neighbourMass / _neighbourMassDensity) || std::isnan(_gradPosDifference[0]) || std::isnan(_gradPosDifference[1]) || std::isnan(_gradPosDifference[2]))
        {
            continue;
        }

        if (_neighbourMassDensity != 0)
        {
            _pressureForceTemp.set(((getParticle(_id).getPressure() + _neighbourPressure)/2) *
                                    (_neighbourMass / _neighbourMassDensity) * 
                                    _gradPosDifference[0], 
                                    ((getParticle(_id).getPressure() + _neighbourPressure)/2) *
                                    (_neighbourMass / _neighbourMassDensity) * 
                                    _gradPosDifference[1],
                                    ((getParticle(_id).getPressure()+ _neighbourPressure)/2) *
                                    (_neighbourMass / _neighbourMassDensity) * 
                                    _gradPosDifference[2]);
        }
        else
        {
            _pressureForceTemp.set(0, 0, 0);
        }
        _pressureForceFinal += _pressureForceTemp;
    }
    getParticle(_id).setPressureForce(-_pressureForceFinal);
}

void System::updateViscosity(const unsigned int &_id)
{
    ngl::Vec3 viscosityFinal = {0,0,0};
    float neighbourMass = 0;
    float neighbourMassDensity = 0;
    ngl::Vec3 neighbourPos = {0,0,0};
    ngl::Vec3 posDifference = {0,0,0};
    ngl::Vec3 neighbourVelocity = {0,0,0};
    ngl::Vec3 velocityDifference = {0,0,0};

    for(const auto &i : getParticle(_id).getNeighbours())
    {
        neighbourMass = getParticle(i).getMass();
        neighbourMassDensity = getParticle(i).getDensity();
        neighbourPos = getParticle(i).getPosition();
        posDifference = getParticle(_id).getPosition() - neighbourPos; 
        neighbourVelocity = getParticle(i).getVelocity();
        velocityDifference = getParticle(i).getVelocity() - neighbourVelocity;
        ngl::Vec3 viscosityTemp(0,0,0);

        if(neighbourMassDensity == 0 || neighbourMass == 0)
        {
            //std::cerr << "Error in updateViscosity: NaN value detected." << std::endl;
            return;
        }
        
        viscosityTemp.set(velocityDifference[0] * (neighbourMass / neighbourMassDensity) *
                          kernelLaplacian(posDifference, 2),
                          velocityDifference[1] * (neighbourMass / neighbourMassDensity) *
                          kernelLaplacian(posDifference, 2),
                          velocityDifference[2] * (neighbourMass / neighbourMassDensity) *
                          kernelLaplacian(posDifference, 2));
        
        if (std::isnan(viscosityTemp.m_x) || std::isnan(viscosityTemp.m_y) || std::isnan(viscosityTemp.m_z))
        {
            //std::cerr << "Error in updateViscosity: NaN value detected." << std::endl;
            return;
        }

        viscosityFinal += viscosityTemp;
    }

    viscosityFinal *= globs._mvisc;

    if (std::isnan(viscosityFinal.m_x) || std::isnan(viscosityFinal.m_y) || std::isnan(viscosityFinal.m_z))
    {
        //std::cerr << "Error in updateViscosity: NaN value detected." << std::endl;
        return;
    }

    getParticle(_id).setViscosity(viscosityFinal);
}

void System::updateGravityForce(const unsigned int &_id)
{
    ngl::Vec3 _gravityForce = getParticle(_id).getMass() * getParticle(_id).getGravityForce();
    getParticle(_id).setGravityForce(_gravityForce);
}

void System::updateBuoyancyForce(const unsigned int &_id)
{
    ngl::Vec3 _buoyancy = {0,0,0};
    _buoyancy.set(globs._buoyancyConst*(getParticle(_id).getDensity() - globs._rd) * getParticle(_id).getGravityForce()[0],
                  globs._buoyancyConst*(getParticle(_id).getDensity() - globs._rd) * getParticle(_id).getGravityForce()[1],
                  globs._buoyancyConst*(getParticle(_id).getDensity() - globs._rd) * getParticle(_id).getGravityForce()[2]);
    getParticle(_id).setBuoyancyForce(_buoyancy);
}

void System::updateNormalField(const unsigned int &_id)
{
    static ngl::Vec3 _normalField = {0, 0, 0};
    float _neighbourMass = 0;
    float _neighbourMassDensity = 0;
    ngl::Vec3 _posDifference = {0,0,0};
    ngl::Vec3 _gradPoly6 = {0,0,0};
    ngl::Vec3 _normalFieldTemp = {0,0,0};

    for(const auto &i: getParticle(_id).getNeighbours())
    {
        float _neighbourMass = getParticle(i).getMass();
        float _neighbourMassDensity = getParticle(i).getDensity();
        if (_neighbourMassDensity <= 0 || std::isnan(_neighbourMassDensity))
        {
            _normalField.set(0,0,0);
            break;
        }
        _posDifference = getParticle(_id).getPosition() - getParticle(i).getPosition();
        if (_posDifference.length() == 0 || std::isnan(_posDifference.length()))
        {
            continue;
        }
        _gradPoly6 = kernelGradient(_posDifference, 0);
        _normalFieldTemp.set(_neighbourMass/ _neighbourMassDensity * _gradPoly6[0],
                             _neighbourMass/ _neighbourMassDensity * _gradPoly6[1],
                             _neighbourMass/ _neighbourMassDensity * _gradPoly6[2]);
        _normalField += _normalFieldTemp;
    }
    getParticle(_id).setNormalField(_normalField); 
}

void System::updateLaplNormalField(const unsigned int &_id)
{
    float _normalLaplField = 0;
    float _neighbourMass = 0;
    float _neighbourMassDensity = 0;
    ngl::Vec3 _neighbourPos = {0,0,0};
    ngl::Vec3 _posDifference = {0,0,0};
    float _normalLaplFieldTemp = 0;

    for(const auto &i: getParticle(_id).getNeighbours())
    {
        _neighbourMass = getParticle(i).getMass();
        _neighbourMassDensity = getParticle(i).getDensity();
        _neighbourPos = getParticle(i).getPosition();
        _posDifference = getParticle(_id).getPosition() - _neighbourPos;
        _normalLaplFieldTemp =  (_neighbourMass/ _neighbourMassDensity) * kernelLaplacian(_posDifference, 0);
        if (std::isnan(_normalLaplFieldTemp)) {
            //std::cerr << "Error: Laplacian of normal field is NaN at particle " << _id << std::endl;
            continue;
        }
        _normalLaplField += _normalLaplFieldTemp;
    }
    getParticle(_id).setLaplNormalField(_normalLaplField);
}

void System::updateSurfaceCurvature(const unsigned int &_id)
{
    updateLaplNormalField(_id);
    updateNormalField(_id);
    ngl::Vec3 normalField = {0,0,0};
    float _surfaceCurvature = 0;

    normalField = getParticle(_id).getNormalField();
    float normalFieldLength = normalField.length();
    if (std::isnan(normalField.m_x) || std::isnan(normalField.m_y) || std::isnan(normalField.m_z) || std::isnan(normalFieldLength) || std::isinf(normalFieldLength) || normalFieldLength == 0)
    {
        //std::cerr << "Error: Normal field is NaN or invalid at particle " << _id << std::endl;
        return;
    }
    _surfaceCurvature = -1 * getParticle(_id).getLaplNormalField() / normalFieldLength;
    if (std::isnan(_surfaceCurvature) || std::isinf(_surfaceCurvature))
    {
        std::cerr << "Error: Surface curvature is NaN or infinite at particle " << _id << std::endl;
        return;
    }
    getParticle(_id).setSurfaceCurvature(_surfaceCurvature);
}


void System::updateSurfaceTension(const unsigned int &_id)
{
    updateSurfaceCurvature(_id);
    ngl::Vec3 _surfaceTension = {0,0,0};
    if((getParticle(_id).getNormalField().length())>=globs._lsurf)
    {
        _surfaceTension.set(getParticle(_id).getSurfaceCurvature()*getParticle(_id).getNormalField()[0],
                            getParticle(_id).getSurfaceCurvature()*getParticle(_id).getNormalField()[1],
                            getParticle(_id).getSurfaceCurvature()*getParticle(_id).getNormalField()[2]);
        getParticle(_id).setSurfaceTension(_surfaceTension);
    }
}

void System::TankCollision(const unsigned int &_id)
{
    ngl::Vec3 _vel = {0,0,0};
    _vel = getParticle(_id).getVelocity();
    if(getParticle(_id).getPosition()[0]>=(m_tank.getDim()[0]+m_tank.getPos()[0]) ||
       getParticle(_id).getPosition()[0]<=(-1*m_tank.getDim()[0]+m_tank.getPos()[0]) )
    {
        _vel[0] *= -1;
    }
    else if(getParticle(_id).getPosition()[1]>=(m_tank.getDim()[1]+m_tank.getPos()[1]) ||
            getParticle(_id).getPosition()[1]<=(-1*m_tank.getDim()[1]+m_tank.getPos()[1]) )
    {
        _vel[1] *= -1;
    }
    else if(getParticle(_id).getPosition()[2]>=(m_tank.getDim()[2]+m_tank.getPos()[2]) ||
            getParticle(_id).getPosition()[2]<=(-1*m_tank.getDim()[2]+m_tank.getPos()[2]) )
    {
        _vel[2] *= -1;
    }
    getParticle(_id).setVelocity(_vel);
}


//UI CHANGES
void System::setBuoyancy(float buoyancyConst)
{
    std::cout<<"Buoyancy set"<<std::endl;
    globs._buoyancyConst = buoyancyConst;
}

void System::setPressure(float _pressureConst)
{
    std::cout<<"Pressure set"<<std::endl;
    globs._kconst = _pressureConst;
}

void System::setGravity(float _fg)
{
    ngl::Vec3 _gravity= {0, 0, 0};
    for (auto &p: m_particles)
    {
        _gravity.m_y = static_cast<float>(_fg);
        p.setGravityForce(_gravity);
    }
    std::cout<<"Gravity set"<<std::endl;
}


void System::setCellSize(float _h)
{
    globs._h =_h;
    std::cout<<"cell size set"<<std::endl;
}

void System::setParticleSize(int _particleSize)
{
    globs.particleSize = _particleSize;
    std::cout<<"Buoyancy set"<<std::endl;
}


void System::setParticleSep(float _particleSep)
{
    globs.particleSep = _particleSep;
    std::cout<<"particle sep set"<<std::endl;
}


void System::setParticleXDim(float _dimX)
{
    globs._dims.m_x = _dimX;
    std::cout<<"particle dim x set"<<std::endl;
}

void System::setParticleYDim(float _dimY)
{
    globs._dims.m_y = _dimY;
    std::cout<<"particle dim y set"<<std::endl;
}

void System::setParticleZDim(float _dimZ)
{
    globs._dims.m_z = _dimZ;
    std::cout<<"particle dim z set"<<std::endl;
}

void System::setRestDensity(float _rd)
{
    globs._rd = _rd;
    std::cout<<"rest density set"<<std::endl;
}


void System::setViscosity(float _viscosity)
{
    globs._mvisc = _viscosity;
    std::cout<<"viscosity set"<<std::endl;
}

void System::setMass(float _mass)
{
    globs.mass = _mass;
    std::cout<<"mass set"<<std::endl;
}


void System::setRLOS(float _RLOS)
{
    globs.RLOS = _RLOS;
    std::cout<<"RLOS set"<<std::endl;
}


void System::setTotalParticles(int _total_particles)
{
    globs.total_particles = _total_particles;
    std::cout<<"Total Particles set"<<std::endl;
}

void System::setTankDimX(float _tankDimX)
{
    m_tank.m_TankDim.m_x = _tankDimX;
    std::cout<<"tank dim x set"<<std::endl;
}

void System::setTankDimY(float _tankDimY)
{
    m_tank.m_TankDim.m_y = _tankDimY;
    std::cout<<"tank dim y set"<<std::endl;
}

void System::setTankDimZ(float _tankDimZ)
{
    m_tank.m_TankDim.m_z = _tankDimZ;
    std::cout<<"tank dim z set"<<std::endl;
}

void System::setTankPosX(float _tankPosX)
{
    m_tank.m_TankPos.m_x = _tankPosX;
    std::cout<<"tank pos x set"<<std::endl;
}

void System::setTankPosY(float _tankPosY)
{
    m_tank.m_TankPos.m_y = _tankPosY;
    std::cout<<"tank pos y set"<<std::endl;
}

void System::setTankPosZ(float _tankPosZ)
{
    m_tank.m_TankPos.m_z = _tankPosZ;
    std::cout<<"tank pos z set"<<std::endl;
}

void System::setDeltaTime(float _dt)
{
    globs.dt = _dt;
    std::cout<<"delta time set"<<std::endl;
}

void System::clearAllParticles()
{
    m_clearAllParticles = true;
    m_particles.clear();
}

void System::ResetSimulation()
{
    m_particles.clear();
    particleReset(globs.total_particles);
    m_ResetSimulation = true;
}

void System::startSimulation()
{
    if(m_particles.size() != 0)
    {
        m_particles.clear();
    }
    particleReset(globs.total_particles);
    m_StartSimulation = true;
}

