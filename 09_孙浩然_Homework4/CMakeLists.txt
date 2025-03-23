# Get all directories
file(GLOB entries RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/* )
set(dirs "")
foreach (entry IN LISTS entries)
    if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${entry} )
        list(APPEND dirs ${CMAKE_CURRENT_SOURCE_DIR}/${entry})
    endif ()
endforeach ()
message(STATUS "[Submission Nodes] " ${dirs})
# Get all targets
set(inner_targets "")
foreach (dir IN LISTS dirs)
    message(STATUS "[Add directory] " ${dir})
    add_subdirectory(${dir})
    get_property(current_targets DIRECTORY ${dir}/nodes PROPERTY BUILDSYSTEM_TARGETS)
    message(STATUS "[Sub targets] " ${current_targets})
    foreach (tgt IN LISTS current_targets)
        message(STATUS "[Pushback target]  " ${tgt})
        list(APPEND inner_targets ${tgt})
    endforeach ()
endforeach ()

message(STATUS "[All sub targets] " ${inner_targets})

get_filename_component(submission_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_custom_target(${submission_name}_nodes
    ALL
        COMMENT "Building submissions..."
)
add_dependencies(${submission_name}_nodes ${inner_targets})

