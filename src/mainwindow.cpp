#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_scene = new NGLScene(this);
    ui->m_mainWindowGridLayout->addWidget(m_scene, 0, 0, 9,2);

    connect(ui->m_massSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setMassFunc(double)));
    connect(ui->m_buoyancySpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setBuoyancyFunc(double)));
    connect(ui->m_RLOSSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setRLOSFunc(double)));
    connect(ui->m_CellSizeSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setCellSizeFunc(double)));
    connect(ui->m_ParticleSepSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setParticleSepFunc(double)));
    connect(ui->m_ParticleSizeSpin, SIGNAL(valueChanged(int)), m_scene, SLOT(setParticleSizeFunc(int)));
    connect(ui->m_RestDensitySpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setRestDensityFunc(double)));
    connect(ui->m_pressureSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setPressureFunc(double)));
    connect(ui->m_GravitySpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setGravityFunc(double)));
    connect(ui->m_NoParticlesSpin, SIGNAL(valueChanged(int)), m_scene, SLOT(setTotalParticlesFunc(int)));
    connect(ui->m_deltaSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setDeltaTimeFunc(double)));

    connect(ui->m_PosXSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setParticleXDimFunc(double)));
    connect(ui->m_PosYSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setParticleYDimFunc(double)));
    connect(ui->m_PosZSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setParticleZDimFunc(double)));

    connect(ui->m_TankDimXSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankDimXFunc(double)));
    connect(ui->m_TankDimYSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankDimYFunc(double)));
    connect(ui->m_TankDimZSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankDimZFunc(double)));

    connect(ui->m_TankPosXSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankPosXFunc(double)));
    connect(ui->m_TankPosYSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankDimYFunc(double)));
    connect(ui->m_TankPosZSpin, SIGNAL(valueChanged(double)), m_scene, SLOT(setTankDimZFunc(double)));

    connect(ui->m_StartSimButton, SIGNAL(clicked(bool)), m_scene, SLOT(StartSimulationFunc()));
    connect(ui->m_stopSimulation, SIGNAL(clicked(bool)), m_scene, SLOT(StopSimulationFunc()));
    connect(ui->m_ResetSimButton, SIGNAL(clicked(bool)), m_scene, SLOT(ResetSimulationFunc()));
    connect(ui->m_ClearSimButton, SIGNAL(clicked(bool)), m_scene, SLOT(clearAllParticlesFunc()));
}

MainWindow::~MainWindow()
{
    delete m_scene;
    delete ui;
}

