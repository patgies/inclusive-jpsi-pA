
add_library(
	amplitude
	STATIC
	amplitudelib.cpp
	datafile.cpp
	tools.cpp
	interpolation.cpp
)

add_executable(
        dipole
        integrand.cpp
        main.cpp
)	

target_link_libraries(
	dipole
	PRIVATE
	cuba
	amplitude
	GSL::gsl
	GSL::gslcblas
)
