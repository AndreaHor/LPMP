add_library(ldp_directed_graph ldp_directed_graph.cxx)
target_link_libraries(ldp_directed_graph LPMP pybind11::module)

add_library(ldp_paths_extractor ldp_paths_extractor.cxx)
target_link_libraries(ldp_paths_extractor ldp_directed_graph LPMP pybind11::module)

add_library(ldp_interval_connection ldp_interval_connection.cxx)
target_link_libraries(ldp_interval_connection ldp_directed_graph ldp_paths_extractor LPMP pybind11::module)

add_library(ldp_instance ldp_instance.cxx )

target_link_libraries(ldp_instance LPMP  pybind11::module ldp_directed_graph ldp_paths_extractor ldp_interval_connection)

add_library(ldp_two_layer_graph ldp_two_layer_graph.cxx)
target_link_libraries(ldp_two_layer_graph LPMP)

add_library(ldp_cut_factor ldp_cut_factor.cxx)
target_link_libraries(ldp_cut_factor LPMP ldp_two_layer_graph)

add_library(ldp_path_factor ldp_path_factor.cxx)
target_link_libraries(ldp_path_factor LPMP pybind11::module)


add_library(ldp_batch_process ldp_batch_process.cxx)
target_link_libraries(ldp_batch_process LPMP pybind11::module)
 
add_library(ldp_fw ldp_fw.cxx)
target_link_libraries(ldp_fw LPMP pybind11::module)


#add_library(lifted_disjoint_paths_input lifted_disjoint_paths_input.cpp)
#target_link_libraries(lifted_disjoint_paths_input LPMP ldp_instance pybind11::module)

add_executable(lifted_disjoint_paths_text_input  lifted_disjoint_paths_text_input.cpp)
target_link_libraries(lifted_disjoint_paths_text_input ldp_instance ldp_paths_extractor ldp_interval_connection ldp_cut_factor ldp_path_factor ldp_directed_graph ldp_batch_process LPMP pybind11::module)


add_executable(lifted_disjoint_paths_fw ldp_fw_solve.cpp)
target_link_libraries(lifted_disjoint_paths_fw ldp_instance ldp_interval_connection ldp_directed_graph ldp_batch_process ldp_fw LPMP FW-MAP pybind11::module)

add_executable(output_mc_from_ldp output_mc_from_ldp.cxx)
target_link_libraries(output_mc_from_ldp ldp_directed_graph LPMP pybind11::module)


pybind11_add_module(ldpMessagePassingPy ldp_python.cxx)

target_link_libraries(ldpMessagePassingPy PRIVATE ldp_instance ldp_interval_connection ldp_paths_extractor ldp_cut_factor ldp_path_factor ldp_directed_graph  ldp_batch_process LPMP pybind11::module)

configure_file(solveFromFiles.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromFiles.py)
#configure_file(solveFromFilesForBatch.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromFilesForBatch.py)
configure_file(solveFromVectors.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromVectors.py)
configure_file(solveFromVectorsTwoGraphs.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromVectorsTwoGraphs.py)
configure_file(solveInBatches.py ${CMAKE_CURRENT_BINARY_DIR}/solveInBatches.py)
#configure_file(solveFromExamples.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromExamples.py)



