set(deps
		stage nodes_system usd geometry usdShade Eigen3::Eigen autodiff igl::core igl_restricted::triangle
)
# If the student add a custom utility lib
message(STATUS "[Utils status: ${student_name}_${util_lib_target_name}_AVAILABLE] " ${${student_name}_${util_lib_target_name}_AVAILABLE})
if (${student_name}_${util_lib_target_name}_AVAILABLE)
	message(STATUS "[${student_name}_${util_lib_target_name}_AVAILABLE] " Available!)
	list(APPEND deps ${student_name}_${util_lib_target_name})
endif ()

add_nodes_with_prefix(
	TARGET_NAME nodes_hw_submissions
	CUS_PREFIX ${student_name}
	DEP_LIBS ${deps}
	COMPILE_DEFS NOMINMAX CGHW_STUDENT_NAME=${student_name}
		EXTRA_INCLUDE_DIRS ${editor_dir} ${editor_dir}/geometry_nodes
)
