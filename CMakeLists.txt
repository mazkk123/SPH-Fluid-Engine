cmake_minimum_required(VERSION 3.5)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})

   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})

endif()


project(FluidGUIBUILD VERSION 0.1 LANGUAGES CXX)


set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)

set(CMAKE_AUTOUIC_SEARCH_PATHS ${PROJECT_SOURCE_DIR}/ui)

set(CMAKE_AUTOMOC ON)

set(CMAKE_AUTORCC ON)


set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(NGL CONFIG REQUIRED)



include_directories(${PROJECT_SOURCE_DIR}/include)



find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(FluidGUI)


target_sources(FluidGUI PRIVATE


    src/main.cpp

    src/mainwindow.cpp

    src/NGLScene.cpp

    src/Tank.cpp

    src/Globals.cpp

    src/System.cpp

    src/Particle.cpp

    src/NGLSceneMouseControls.cpp

    include/mainwindow.h

    include/NGLScene.h

    include/WindowParams.h

    include/Tank.h

    include/Globals.h

    include/System.h

    include/Particle.h

    ui/mainwindow.ui


)



target_link_libraries(FluidGUI PRIVATE NGL Qt5::Widgets)

if (Qt6_FOUND)

    target_link_libraries(FluidGUI PRIVATE Qt::OpenGLWidgets)

endif()


add_custom_target(${TargetName}CopyShadersAndfonts ALL

    COMMAND ${CMAKE_COMMAND} -E copy_directory

    ${CMAKE_CURRENT_SOURCE_DIR}/shaders

    $<TARGET_FILE_DIR:FluidGUI>/shaders

)
