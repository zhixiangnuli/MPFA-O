cmake_minimum_required(VERSION 3.16)
project(Tensor_MPFA)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 23)
set(MKL_LINK "static")
set(MKL_INTERFACE "lp64")
find_package(MKL)

add_executable(Tensor_MPFA TEST.cpp GMRES.cpp GMRES.h Checked.h tran.h roots.h utility.h)
target_include_directories(Tensor_MPFA PUBLIC thirdparty\\eigen3 thirdparty\\mdspan\\include\\experimental )
target_include_directories(Tensor_MPFA SYSTEM PUBLIC 
#	$<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
	"${MKL_H}"
)
target_link_libraries(Tensor_MPFA PRIVATE
	MKL::MKL
)