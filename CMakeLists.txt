#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.16)

#Name your project here
project(3D_VDB_STENCIL LANGUAGES CUDA CXX)

#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_CXX_STANDARD 20)

# fix warnings with BOOST and CMake
#set(Boost_NO_WARN_NEW_VERSIONS 1)
#find_package(Boost REQUIRED COMPONENTS program_options)


find_package(OpenVDB CONFIG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(TIRA REQUIRED)
#find_package(OpenCV REQUIRED)


#build the executable in the binary directory on MS Visual Studio
if ( MSVC )
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
	SET( LIBRARY_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( LIBRARY_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif ( MSVC )



#set the include directories
include_directories(
			${CMAKE_CURRENT_BINARY_DIR}
			${CMAKE_CURRENT_SOURCE_DIR}
			${TIRA_INCLUDE_DIRS}
                         
)
link_directories(${CUDA_BIN_DIRS})
				
#create an executable
add_executable(3D_VDB_STENCIL
					main.cpp
)


target_link_libraries(3D_VDB_STENCIL
						PRIVATE OpenVDB::openvdb
)