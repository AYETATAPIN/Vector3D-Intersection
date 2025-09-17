set(core_SRCS ledas_task.cpp)

add_library(core ${core_SRCS})
add_executable(exe ledas_task.cpp)
target_link_libraries(exe core)