#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <ngl/Vec3.h>
#include <ngl/AbstractVAO.h>
#include <list>
#include <unordered_map>
#include "Globals.h"
#include "Particle.h"

class System
{
    public:
        System(ngl::Vec3 _pos={0.0f,0.0f,0.0f}, unsigned int frame = 1);
        ~System();

        //QT METHODS
        
        void setGravity(float _fg);

        void setPressure(float _pressureConst);

        void setBuoyancy(float buoyancyConst);

        void setCellSize(float _h);

        void setParticleSize(int _particleSize);

        void setParticleSep(float _particleSep);

        void setParticleXDim(float _dimX);
        void setParticleYDim(float _dimY);
        void setParticleZDim(float _dimZ);

        void setRestDensity(float _rd);

        void setViscosity(float _viscosity);

        void setMass(float _mass);

        void setRLOS(float _RLOS);

        void setTotalParticles(int _total_particles);

        void setTankDimX(float _tankDimX);
        void setTankDimY(float _tankDimY);
        void setTankDimZ(float _tankDimZ);

        void setTankPosX(float _tankPosX);
        void setTankPosY(float _tankPosY);
        void setTankPosZ(float _tankPosZ);

        void setDeltaTime(float _dt);

        void clearAllParticles();
        void ResetSimulation();
        void startSimulation();


        //HASH TABLE ASSIGNS

        void update();
        void draw();

        void addToHashTable();
        void updateHashTable(Particle &_p);
        void systemRenderFrames();
        bool isPrime(unsigned int &_numToCheck);
        ngl::Vec3 x_hat(const ngl::Vec3 &m_pos);
        unsigned int primeFunction(unsigned int _n);
        Particle &getParticle(const unsigned int &_i);

        void particleReset(unsigned int &_numParticles);
        void resetParticleToTank();
        void resetParticleToRandom(Particle &_p);
        void resetParticleToSphere();

	    ngl::Vec3 randomVectorOnSphere(float _radius = 1.0f);

        unsigned int XOR(const unsigned int &_a, const unsigned int &_b);
        std::unordered_map<unsigned int, std::vector<unsigned int>>getHashTable();

        bool m_export;
        std::unordered_map<unsigned int, std::vector<unsigned int>>m_hashMap;
        std::vector<Particle>m_particles;
        std::vector<int>m_neighbours;

        //PARTICLE CALCULATIONS 

        float kernel(const ngl::Vec3 &_r, const int &_type);
        ngl::Vec3 kernelGradient(const ngl::Vec3 &_r, const int &_type);
        float kernelLaplacian(const ngl::Vec3 &_r, const int &_type);

        void findNeighbours(const unsigned int &_id);
        void updatePosition(const unsigned int &_id);
        void updateMass(unsigned int &_id);
        void updateVelocity(const unsigned int &_id);
        void updateXSPHVelocity(const unsigned int &_id);
        void updateAcceleration(const unsigned int &_id);
        void updateDensity(const unsigned int &_id);
        void updatePressure(const unsigned int &_id);
        void updatePressureForce(const unsigned int &_id);
        void updateViscosity(const unsigned int &_id);
        void updateGravityForce(const unsigned int &_id);
        void updateBuoyancyForce(const unsigned int &_id);
        void updateSumForces(const unsigned int &_id);
        void updateNormalField(const unsigned int &_id);
        void updateLaplNormalField(const unsigned int &_id);
        void updateSurfaceCurvature(const unsigned int &_id);
        void updateSurfaceTension(const unsigned int &_id);
        void TankCollision(const unsigned int &_id);

        void updateCurrFrame();
        void updateTotalFrames();

        Globals globs;
        Tank m_tank;
        unsigned int m_numParticles;
        ngl::Vec3 m_pos;
        std::vector<unsigned int>ids={};
        std::vector<int>unwanted_data={};
        unsigned int m_frame = 1;
        unsigned int m_currFrame = 0;

        bool m_clearAllParticles;
        bool m_ResetSimulation;
        bool m_StartSimulation;

        std::unique_ptr<ngl::AbstractVAO> m_vao;
};

#endif
