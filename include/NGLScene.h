#ifndef NGLSCENE_H_
#define NGLSCENE_H_
#include <ngl/Vec3.h>
#include "WindowParams.h"
#include <memory>
#include "System.h"
#include "Particle.h"
#include "Globals.h"
#include "Tank.h"
// this must be included after NGL includes else we get a clash with gl libs
#include <QOpenGLWidget>
//----------------------------------------------------------------------------------------------------------------------
/// @file NGLScene.h
/// @brief this class inherits from the Qt OpenGLWindow and allows us to use NGL to draw OpenGL
/// @author Jonathan Macey
/// @version 1.0
/// @date 10/9/13
/// Revision History :
/// This is an initial version used for the new NGL6 / Qt 5 demos
/// @class NGLScene
/// @brief our main glwindow widget for NGL applications all drawing elements are
/// put in this file
//----------------------------------------------------------------------------------------------------------------------

class NGLScene : public QOpenGLWidget
{
   Q_OBJECT
  public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor for our NGL drawing class
    /// @param [in] parent the parent window to the class
    //----------------------------------------------------------------------------------------------------------------------
    NGLScene(QWidget *_parent);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor must close down ngl and release OpenGL resources
    //----------------------------------------------------------------------------------------------------------------------
    ~NGLScene() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the initialize class is called once when the window is created and we have a valid GL context
    /// use this to setup any default GL stuff
    //----------------------------------------------------------------------------------------------------------------------
    void initializeGL() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we want to draw the scene
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we resize the window
    //----------------------------------------------------------------------------------------------------------------------
    void resizeGL(int _w, int _h) override;

public slots:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------
    void setBuoyancyFunc(double _buoyancyConst);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system using a slider
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setCellSizeFunc(double _h);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of cell size in the system using a slider
    /// @param [in] float the cell size value to be set with

    void setParticleSizeFunc(int _particleSize);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of individual particle sizes' in the system from a slider
    /// @param [in] float the particle size value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setParticleSepFunc(double _particleSep);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setParticleXDimFunc(double _dimX);
    void setParticleYDimFunc(double _dimY);
    void setParticleZDimFunc(double _dimZ);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------
    void setRestDensityFunc(double _rd);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setPressureFunc(double _fp);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setViscosityFunc(double _viscosity);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setMassFunc(double _mass);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setRLOSFunc(double _RLOS);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setGravityFunc(double _fg);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setTotalParticlesFunc(int _total_particles);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------

    void setTankDimXFunc(double _tankDimX);
    void setTankDimYFunc(double _tankDimY);
    void setTankDimZFunc(double _tankDimZ);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------
    void setTankPosXFunc(double _tankPosX);
    void setTankPosYFunc(double _tankPosY);
    void setTankPosZFunc(double _tankPosZ);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the value of buoyancy force in the system
    /// @param [in] float the buoyancy value to be set with
    //----------------------------------------------------------------------------------------------------------------------
    void setDeltaTimeFunc(double _dt);

    void clearAllParticlesFunc();
    void ResetSimulationFunc();
    void StartSimulationFunc();
    void StopSimulationFunc();

private:

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Qt Event called when a key is pressed
    /// @param [in] _event the Qt event to query for size etc
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *_event) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called every time a mouse is moved
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseMoveEvent (QMouseEvent * _event ) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is pressed
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mousePressEvent ( QMouseEvent *_event) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is released
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseReleaseEvent ( QMouseEvent *_event ) override;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse wheel is moved
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void wheelEvent( QWheelEvent *_event) override;
    /// @brief windows parameters for mouse control etc.

    void timerEvent(QTimerEvent *) override;

    WinParams m_win;
    /// position for our model
    ngl::Vec3 m_modelPos;
    ngl::Mat4 m_view;
    ngl::Mat4 m_project;
    ngl::Mat4 m_globalMouseTX;
    int m_timer = 0;
    std::unique_ptr<System> m_sys;
    bool m_stopSim = false;
    std::chrono::high_resolution_clock::time_point m_previousTime;

};



#endif
