REMARK: Using scotch reduce the execution time of mmg3d (see https://gforge.inria.fr/frs/?group_id=248 to download it, mmg3d work with Scotch up to 5.1.12).
 If you have already installed Scotch and want to use it, please set
 the CMake variable or environment variable SCOTCH_DIR to your scotch
 directory ( export SCOTCH_DIR=~/scotch_5.0 ).

DIRECTORIES
  sources     = source files of project Mmg3d4

BUILDING AND INSTALLING MMG3D4:
  $cd mmg3d4
  $mkdir build
  $cd build

  FIRST STEP: CONFIGURATION
      $ccmake ..
      # - you must configure (c) then if needed you can toogle to advanced mode (t).
      #   To end, generate the makefile with (g)
    #or:
      $cmake ..

      # MAIN OPTIONS:
      # - CMAKE_INSTALL_PREFIX allow to configure the location of the runtimes (default /usr/local).
      # - You can (dis)able the compilation of shared (LIBMMG3D4_SHARED) and static (LIBMMG3D4_STATIC) mmg3d4 libraries and the compilation of examples of use of the libraries (TEST_LIBMMG3D4 variable).
      # - You can (dis)able the use of scotch (USE_SCOTCH). Work only with 5.1.12 release
      # - variable CMAKE_BUILD_TYPE can be set to Debug (for -g flag), Release (for -O3 flag)

  SECOND STEP: COMPILATION
      $make

  THIRD STEP: INSTALLATION
      $make install
      # Header files are installed in ${CMAKE_INSTALL_PREFIX}/include
      # Libraries in ${CMAKE_INSTALL_PREFIX}/lib
      # Runtimes in ${CMAKE_INSTALL_PREFIX}/bin

OUTPUTS:
  mmg3d4       in normal mode
  mmg3d4_debug in Debug mode
  mmg3d4_O3    in Release mode
  mmg3d4_Os    in MinSizeRel mode
  if activated :
     libmmg3d4.a
     libmmg3d4.so
