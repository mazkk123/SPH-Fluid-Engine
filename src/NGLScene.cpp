#include <QMouseEvent>
#include <QGuiApplication>

#include "NGLScene.h"
#include <ngl/NGLInit.h>
#include <ngl/Util.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOFactory.h>
#include <ngl/VAOPrimitives.h>
#include <iostream>

NGLScene::NGLScene(QWidget *_parent) : QOpenGLWidget(_parent)
{
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  //setTitle("Fluid System");
}

NGLScene::~NGLScene()
{
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
}

void NGLScene::resizeGL(int _w , int _h)
{
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
  m_project = ngl::perspective(45.0f, static_cast<float>(_w)/_h, 0.1f,600.0f );
}

void NGLScene::initializeGL()
{
  // we must call that first before any other GL commands to load and link the
  // gl commands from the lib, if that is not done program will crash
  ngl::NGLInit::initialize();
  m_view = ngl::lookAt({0.4f, 0.4f, 0.4f}, {0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
  glClearColor(0.2f, 0.2f, 0.2f, 1.0f);			   // Grey Background
  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);
  // enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);

  //ngl::VAOPrimitives::createTrianglePlane("ground",4,4,1,1,ngl::Vec3::up());

  ngl::ShaderLib::use(ngl::nglColourShader);
  ngl::ShaderLib::setUniform("Colour",0.6f,0.6f,0.6f,1.0f);

  ngl::ShaderLib::loadShader("ParticleShader","shaders/ParticleVertex.glsl", "shaders/ParticleFragment.glsl");

  m_sys = std::make_unique<System>();
  startTimer(0);
  m_previousTime = std::chrono::high_resolution_clock::now();
  glPointSize(m_sys->globs.particleSize);

}

void NGLScene::paintGL()
{
  // clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_win.width,m_win.height);

  auto xrot=ngl::Mat4::rotateX(m_win.spinXFace);
  auto yrot= ngl::Mat4::rotateY(m_win.spinYFace);
  m_globalMouseTX = xrot * yrot;
  m_globalMouseTX.m_m[3][0] = m_modelPos.m_x;
  m_globalMouseTX.m_m[3][1] = m_modelPos.m_y;
  m_globalMouseTX.m_m[3][2] = m_modelPos.m_z;

  ngl::ShaderLib::use(ngl::nglColourShader);
  ngl::ShaderLib::setUniform("MVP",m_project*m_view * m_globalMouseTX);
  ngl::ShaderLib::use("ParticleShader");
  ngl::ShaderLib::setUniform("MVP",m_project*m_view * m_globalMouseTX);

  if(m_sys->m_StartSimulation || m_sys->m_ResetSimulation )
  {
        m_sys->draw();
        m_timer += 1;
        std::cout<<"Frame "<<m_timer<<std::endl;
  }
  else
  {
        m_sys->m_ResetSimulation = false;
        m_sys->m_StartSimulation = false;
  }
}


//----------------------------------------------------------------------------------------------------------------------

void NGLScene::keyPressEvent(QKeyEvent *_event)
{
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
  {
  // escape key to quite
  case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
  case Qt::Key_Space :
      m_win.spinXFace=0;
      m_win.spinYFace=0;
      m_modelPos.set(ngl::Vec3::zero());

  break;
  default : break;
  }
  // finally update the GLWindow and re-draw

    update();
}

void NGLScene::timerEvent(QTimerEvent *)
{

  auto currentTime = std::chrono::high_resolution_clock::now();
  auto delta = std::chrono::duration<float,std::chrono::seconds::period>(currentTime-m_previousTime).count();
  if(m_sys->m_StartSimulation || m_sys->m_ResetSimulation)
  {
        if(!(m_stopSim))
        {
        m_sys->update();
        }
  }
  update();
  m_previousTime=currentTime;
}


//UI WIDGET CALLS
void NGLScene::setBuoyancyFunc(double _buoyancyConst)
{
    m_sys->setBuoyancy(_buoyancyConst);
}

void NGLScene::setPressureFunc(double _fp)
{
    m_sys->setPressure(_fp);
}


void NGLScene::setRLOSFunc(double _RLOS)
{
    m_sys->setRLOS(_RLOS);
}


void NGLScene::setGravityFunc(double _fg)
{
    m_sys->setGravity(_fg);
}


void NGLScene::setCellSizeFunc(double _h)
{
    m_sys->setCellSize(_h);
}


void NGLScene::setParticleSizeFunc(int _particleSize)
{
    m_sys->setParticleSize(_particleSize);
}


void NGLScene::setParticleSepFunc(double _particleSep)
{
    m_sys->setParticleSep(_particleSep);
}

void NGLScene::setParticleXDimFunc(double _dimX)
{
    m_sys->setParticleXDim(_dimX);
}

void NGLScene::setParticleYDimFunc(double _dimY)
{
    m_sys->setParticleYDim(_dimY);
}

void NGLScene::setParticleZDimFunc(double _dimZ)
{
    m_sys->setParticleZDim(_dimZ);
}

void NGLScene::setRestDensityFunc(double _rd)
{
    m_sys->setRestDensity(_rd);
}

void NGLScene::setViscosityFunc(double _viscosity)
{
    m_sys->setViscosity(_viscosity);
}

void NGLScene::setMassFunc(double _mass)
{
    m_sys->setMass(_mass);
}

void NGLScene::setTotalParticlesFunc(int _total_particles)
{
    m_sys->setTotalParticles(_total_particles);
}



void NGLScene::setTankDimXFunc(double _tankDimX)
{
    m_sys->setTankDimX(_tankDimX);
}

void NGLScene::setTankDimYFunc(double _tankDimY)
{
    m_sys->setTankDimX(_tankDimY);
}

void NGLScene::setTankDimZFunc(double _tankDimZ)
{
    m_sys->setTankDimY(_tankDimZ);
}

void NGLScene::setTankPosXFunc(double _tankPosX)
{
    m_sys->setTankDimZ(_tankPosX);
}

void NGLScene::setTankPosYFunc(double _tankPosY)
{
    m_sys->setTankPosY(_tankPosY);
}

void NGLScene::setTankPosZFunc(double _tankPosZ)
{
    m_sys->setTankPosZ(_tankPosZ);
}


void NGLScene::setDeltaTimeFunc(double _dt)
{
    m_sys->setDeltaTime(_dt);
}


void NGLScene::clearAllParticlesFunc()
{
    m_sys->clearAllParticles();
}

void NGLScene::ResetSimulationFunc()
{
    m_stopSim = false;
    m_timer = 0;
    m_sys->ResetSimulation();
}

void NGLScene::StartSimulationFunc()
{
    m_stopSim = false;
    m_timer = 0;
    m_sys->startSimulation();
}

void NGLScene::StopSimulationFunc()
{
    m_stopSim = true;
}
