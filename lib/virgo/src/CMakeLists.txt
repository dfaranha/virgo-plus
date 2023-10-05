link_directories(lib)
link_directories("../../../../HElib/build/helib_pack/lib/")
aux_source_directory(. LIB_SRC)
add_library(virgo_lib ${LIB_SRC})
include_directories("../include")
include_directories("../../../../HElib/build/helib_pack/include/")
target_link_libraries(virgo_lib XKCP helib ntl gmp)

