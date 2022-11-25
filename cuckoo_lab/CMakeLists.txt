add_executable(cuckoo_lab main.cpp)

find_package(libcuckoo CONFIG REQUIRED)
target_link_libraries(cuckoo_lab libcuckoo::libcuckoo)