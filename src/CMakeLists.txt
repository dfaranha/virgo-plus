aux_source_directory(. SRC)
include_directories(.)
include_directories(../lib/virgo/include)
add_executable(virgo_plus_run ${SRC})

if(PC_OPTION STREQUAL "virgo")
    target_link_libraries(virgo_plus_run virgo_lib XKCP)
	target_link_libraries(virgo_plus_run virgo_lib gmp)
endif()
