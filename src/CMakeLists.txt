find_package(ArrayFire)
add_executable(pixDPC main.cpp dpc_support.cpp)

# To use Unified backend, do the following.
# Unified backend lets you choose the backend at runtime
target_link_libraries(pixDPC ArrayFire::af)
