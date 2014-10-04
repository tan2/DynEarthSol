CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)
PROJECT (mmg3d4)

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
#SET(EXECUTABLE_OUTPUT_PATH $ENV{HOME}/bin/$ENV{ARCHI}/)
SET(MMG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/sources/)
SET(EXECUTABLE_OUTPUT_PATH ${MMG_SOURCE_DIR}/../bin)

# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
  # using clang
  SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
  # using icc
  SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()


IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "-Wno-char-subscripts ${CMAKE_C_FLAGS}")
  IF(APPLE)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64 ${CMAKE_EXE_LINKER_FLAGS}")

    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
    #  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF()
ENDIF()

FILE(
  GLOB
  heaser_files
  ${MMG_SOURCE_DIR}/*.h
  )
ADD_CUSTOM_COMMAND(OUTPUT .dummy
  COMMAND touch .dummy
  DEPENDS ${header_files}
  )

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${MMG_SOURCE_DIR})
# Source files
FILE(
  GLOB
  source_files
  ${MMG_SOURCE_DIR}/*.c   ${MMG_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM source_files
  ${MMG_SOURCE_DIR}/mmg3d.c ${MMG_SOURCE_DIR}/mmg3dlib.c ${REMOVE_FILE})
FILE(
  GLOB
  main_file
  ${MMG_SOURCE_DIR}/mmg3d.c
  )
FILE(
  GLOB
  lib_file
  ${MMG_SOURCE_DIR}/mmg3dlib.c
  )

IF ((NOT WIN32) AND (NOT WIN64))
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()

###############################################################################
#####
#####         Scotch
#####
###############################################################################
# add SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
INCLUDE(FindScotch.cmake)

IF(SCOTCH_FOUND)
  OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON)
  IF( USE_SCOTCH )
    # Set flags for building test program
    INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIR})

    SET(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARY})
    SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCHERR_LIBRARY})

    SET(CMAKE_C_FLAGS "-DUSE_SCOTCH ${CMAKE_C_FLAGS}")
    SET(SCOTCH_LINK_FLAGS "-lpthread")
    MESSAGE(STATUS "Compilation with scotch")
    SET( LIBRARIES ${SCOTCH_LINK_FLAGS} ${LIBRARIES})
    SET( LIBRARIES ${SCOTCHERR_LIBRARY} ${LIBRARIES})
    SET( LIBRARIES ${SCOTCH_LIBRARY} ${LIBRARIES})
  ENDIF()
ELSE ( )
  MESSAGE ( WARNING "Scotch library not found:
 Using scotch reduce the execution time of mmg3d (see https://gforge.inria.fr/frs/?group_id=248 to download it).
 If you have already installed Scotch and want to use it, please set the CMake variable or environment variable SCOTCH_DIR to your scotch directory.")
  OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" OFF)
ENDIF ( )

###############################################################################
#####
#####         Required Mode
#####
###############################################################################
# preserve required entites?
OPTION ( REQUIRED "Preserve required entities" OFF)
IF(REQUIRED)
  SET(CMAKE_C_FLAGS "-DREQUIRED ${CMAKE_C_FLAGS}")
  MESSAGE("Required quantities preservation enabled")
ENDIF()

###############################################################################
#####
#####         Compile mmg3d4 libraries
#####
###############################################################################
SET(LIBRARY_OUTPUT_PATH ${MMG_SOURCE_DIR}/../lib)

OPTION ( LIBMMG3D4_STATIC "Compile static library" ON)
IF ( LIBMMG3D4_STATIC )
  ADD_LIBRARY(${PROJECT_NAME}_a  STATIC
    .dummy
    ${source_files} ${lib_file} )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_a PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_a ${LIBRARIES})
  INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h DESTINATION include)
  INSTALL(TARGETS ${PROJECT_NAME}_a
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
OPTION ( LIBMMG3D4_SHARED "Compile dynamic library" OFF)
IF ( LIBMMG3D4_SHARED )
  ADD_LIBRARY(${PROJECT_NAME}_so SHARED
    .dummy
    ${source_files} ${lib_file})
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_so PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_so ${LIBRARIES})
  INSTALL(FILES ${MMG_SOURCE_DIR}/lib${PROJECT_NAME}.h DESTINATION include)
  INSTALL(TARGETS ${PROJECT_NAME}_so
   ARCHIVE DESTINATION
   LIBRARY DESTINATION lib)
ENDIF()

###############################################################################
#####
#####         Compile program to test library
#####
###############################################################################
OPTION ( TEST_LIBMMG3D4 "Compile tests of library" OFF)
IF ( TEST_LIBMMG3D4 )
  IF ( LIBMMG3D4_STATIC )
    ADD_EXECUTABLE(libmmg4_example1 ${MMG_SOURCE_DIR}/../libexamples/main.c)
    TARGET_LINK_LIBRARIES(libmmg4_example1 ${PROJECT_NAME}_a)
    INSTALL(TARGETS libmmg4_example1 RUNTIME DESTINATION bin )
  ELSEIF ( LIBMMG3D5_SHARED )
    ADD_EXECUTABLE(libmmg4_example1  ${MMG_SOURCE_DIR}/../libexamples/main.c)
    TARGET_LINK_LIBRARIES(libmmg4_example1 ${PROJECT_NAME}_so)
    INSTALL(TARGETS libmmg4_example1 RUNTIME DESTINATION bin )
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ()
ENDIF()

###############################################################################
#####
#####         Compile executable
#####
###############################################################################
ADD_EXECUTABLE(${PROJECT_NAME}
  .dummy
  ${source_files} ${main_file})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX _debug) # in debug mode we name the executable mmg3d4_debug
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES RELEASE_POSTFIX _O3) # in Release mode we name the executable mmg3d4_O3
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES MINSIZEREL_POSTFIX _Os)

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
OPTION ( BUILD_TESTING "Enable/Disable continuous integration" OFF )
MARK_AS_ADVANCED( BUILD_TESTING )
IF(BUILD_TESTING)
  ENABLE_TESTING()
  INCLUDE(CTest)

  IF ( CMAKE_COMPILER_IS_GNUCC AND NOT APPLE )
    SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -fprofile-arcs -ftest-coverage --coverage ${CMAKE_C_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g -fprofile-arcs -ftest-coverage --coverage ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF ( )

  # Add runtime that we want to test
  GET_TARGET_PROPERTY(EXECUT ${PROJECT_NAME} LOCATION)
  IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET(EXECUT ${EXECUT}_debug)
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    SET(EXECUT ${EXECUT}_O3)
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Minsizerel")
    SET(EXECUT ${EXECUT}_Os)
  ENDIF()
  SET ( LISTEXEC ${EXECUT} )
  IF ( TEST_LIBMMG3D4 )
    IF ( LIBMMG3D4_STATIC OR LIBMMG3D4_SHARED )
      GET_TARGET_PROPERTY(LIBEXEC1 libmmg4_example1 LOCATION)
      SET( LISTEXEC ${LISTEXEC} ${LIBEXEC1} )
    ENDIF()
  ENDIF()

ENDIF(BUILD_TESTING)


###############################################################################
#####
#####         Create Package
#####
###############################################################################
INCLUDE(InstallRequiredSystemLibraries)
SET( CPACK_SOURCE_GENERATOR TGZ )
SET( CPACK_GENERATOR TGZ )

SET(CPACK_PACKAGE_NAME ${PROJECT_NAME})
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG3D: ISOTROPIC TETRAHEDRAL REMESHER")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey")
SET(CPACK_PACKAGE_CONTACT "cecile.dobrzynski@inria.fr")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION "4.0.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_INSTALL_DIRECTORY}-${CMAKE_SYSTEM})
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}")
SET(CPACK_OUTPUT_FILE_PREFIX "../archives")

SET(CPACK_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.svn;.~;.#" ${CPACK_IGNORE_FILES})

SET(CPACK_SOURCE_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.svn;.~;.#" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";lib/;bin;archives"    ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";README_developper.txt"    ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";CMakeFiles" ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";.#$"        ${CPACK_SOURCE_IGNORE_FILES})
SET(CPACK_SOURCE_IGNORE_FILES ";build*" ${CPACK_SOURCE_IGNORE_FILES} )

IF ((WIN32 OR WIN64) AND NOT UNIX)
  # Warning compilation under windows not tested
  #SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/logo.pdf")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\${PROJECT_NAME}.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} MMG3D: ISOTROPIC TETRAHEDRAL REMESHER")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.math.u-bordeaux1.fr/~dobrzyns/.../mmg3d.php‎")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.math.u-bordeaux1.fr/~dobrzyns/.../mmg3d.php‎")
  SET(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})
  SET(CPACK_NSIS_MODIFY_PATH ON) # Choice to add prog dir to the system PATH
  SET(CPACK_PACKAGING_INSTALL_PREFIX  "${PROJECT_NAME}" CACHE PATH "Where you want to install your package")
  MARK_AS_ADVANCED(CPACK_PACKAGING_INSTALL_PREFIX)
ELSEIF (APPLE )
  SET(CPACK_BINARY_OSX11 "ON")
ELSEIF ( UNIX )
ENDIF()

INCLUDE(CPack)
