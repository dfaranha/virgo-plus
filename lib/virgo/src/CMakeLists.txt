link_directories(lib)
aux_source_directory(. LIB_SRC)
add_library(virgo_lib ${LIB_SRC})
include_directories("../include")
target_link_libraries(virgo_lib XKCP)