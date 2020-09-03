macro(artver)
  execute_process(COMMAND ups active COMMAND grep -v "Active UPS products:" COMMAND grep -e "^art" COMMAND awk "{print $2}" OUTPUT_VARIABLE art_ver)
  set(min_art_ver "v2_00_00")
  set(next_art_ver "v2_06_00")
endmacro()

macro(canvasver)
  execute_process(COMMAND ups active COMMAND grep -v "Active UPS products:" COMMAND grep -e "^canvas" COMMAND awk "{print $2}" OUTPUT_VARIABLE canvas_ver)
  set(min_canvas_ver "v1_06_00")
endmacro()

macro(libsForPlugins)


  # decide on which set of libraries to use based on art version
  artver()
  canvasver()
  if (NOT ${art_ver} EQUAL "" )
    if (${art_ver} STRLESS ${min_art_ver})
      # art pre v2.00.00
      list(INSERT service_lib_list 0
        ${art_Framework_Services_Registry}
        ${art_Persistency_Common}
        ${art_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT module_lib_list 0
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT source_lib_list 0
        ${art_Framework_IO_Sources}
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
    elseif (${art_ver} STRLESS ${next_art_ver})
      # art pre v2.06.00
      list(INSERT service_lib_list 0
        ${art_Framework_Services_Registry}
        ${art_Persistency_Common}
        ${art_Utilities}
        ${canvas_Persistency_Common}
        ${canvas_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT module_lib_list 0
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${canvas_Persistency_Common}
        ${canvas_Persistency_Provenance}
        ${canvas_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT source_lib_list 0
        ${art_Framework_IO_Sources}
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${canvas_Persistency_Common}
        ${canvas_Persistency_Provenance}
        ${canvas_Utilities}
        ${fhiclcpp}
        ${cetlib}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
    else()
      # art post v2.06.00
      list(INSERT service_lib_list 0
        ${art_Framework_Services_Registry}
        ${art_Persistency_Common}
        ${art_Utilities}
        ${canvas}
        ${fhiclcpp}
        ${cetlib}
        ${cetlib_except}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT module_lib_list 0
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${canvas}
        ${fhiclcpp}
        ${cetlib}
        ${cetlib_except}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
      list(INSERT source_lib_list 0
        ${art_Framework_IO_Sources}
        ${art_Framework_Core}
        ${art_Framework_Principal}
        ${art_Persistency_Common}
        ${art_Persistency_Provenance}
        ${art_Utilities}
        ${canvas}
        ${fhiclcpp}
        ${cetlib}
        ${cetlib_except}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
    endif()

  elseif(NOT ${canvas_ver} EQ "")
    # What if art is not set up but canvas is
    if (${canvas_ver} STRGREATER ${min_canvas_ver})
      list(INSERT lib_list 0
        ${canvas}
        ${fhiclcpp}
        ${cetlib}
        ${cetlib_except}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
    else()
      list(INSERT lib_list 0
        ${canvas}
        ${fhiclcpp}
        ${cetlib}
        ${root_Core}
        ${boost_filesystem}
        ${boost_system}
        ${MF_MessageLogger})
    endif()
  endif()

endmacro()
