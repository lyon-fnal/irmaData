# The basic pattern for creating shared libraries should look like:
#
#   add_library(MyProducer_module SHARED MyProducer_module.cc)
#   target_link_library(MyProducer_module <list of library dependencies>)
#
# Any of the studio create-* commands will automatically generate
# these lines for the requested plugin.  For those who do not use the
# "create" commands, or who have more complicated use cases, by-hand
# adjustment may be required.

# "irmaData" lib linking
add_library(irmaData_module SHARED irmaData_module.cc)
target_link_libraries(irmaData_module ${module_lib_list})

#message("${module_lib_list}")

# "ratioeast" lib linking
add_library(ratioeast_module SHARED ratioeast_module.cc)
target_link_libraries(ratioeast_module ${module_lib_list}
                                       ${art_Framework_Services_Optional_TFileService_service}
                                       ${art_Framework_Services_Optional}
                                       ${art_Framework_Services_Registry}
                                       ${art_Framework_IO_Root}
                                       ${root_HistFactory} ${root_Hist} ${root_RIO} ${root_MathCore}
                                       )