add_library(Luna SHARED SchemeBindings.h src/SchemeBindings.cpp)
target_link_libraries(Luna PRIVATE LUNA_LIB)

add_custom_target(LunaRunner_GEN)
add_custom_command(TARGET LunaRunner_GEN
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND LUNA_HOME=${PROJECT_SOURCE_DIR}/projects/Core/luna ${PROJECT_SOURCE_DIR}/projects/Core/luna/util/compiler/lunac Scheme.fa Scheme.ja
)

add_dependencies(Luna LunaRunner_GEN)
