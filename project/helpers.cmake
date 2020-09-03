#cmake macro

#macro for grabbing compiler directories
macro(get_comp_paths)
execute_process(COMMAND type -p g++ OUTPUT_VARIABLE gpp_path)
execute_process(COMMAND type -p gcc OUTPUT_VARIABLE gcc_path)
message(STATUS "gpp path: ${gpp_path}")
message(STATUS "gcc path: ${gcc_path}")
endmacro(get_comp_paths)


#map macro
macro(insert_into_map _key _val)
    set("map_${_key}" "${_val}")
endmacro(insert_into_map)


#here is where we find the include directories
macro(find_inc_dirs)

#getting ups product names
    execute_process(COMMAND ups active COMMAND grep -v "Active ups products:" COMMAND cut -d " " -f 1 COMMAND tr '[a-z]' '[A-Z]' OUTPUT_VARIABLE dirs)
#trimming the end
    string(LENGTH ${dirs} string_length)
    string(SUBSTRING ${dirs} 0 ${string_length}-2 dirs_trimmed)
#making the string of dir names into a list
    string(REPLACE "\n" ";" dirs_list ${dirs_trimmed})
#THIS IS THE VARIABLE WE WILL WANT TO USE IN CMakeList
    set(inc_dirs_list)
#looping over product list to grab inc dirs
    foreach(prod ${dirs_list})
        if(prod STREQUAL "")
            continue()
        endif(prod STREQUAL "")
        list(APPEND inc_dirs_list $ENV{${prod}_INC})
    endforeach(prod)
#message(STATUS "${inc_dirs_list}")

endmacro()

#========================

#here is where we find the lib directories
macro(find_lib_dirs)

    execute_process(COMMAND ups active COMMAND grep -v "Active UPS products:" COMMAND cut -d " " -f 1 COMMAND tr '[a-z]' '[A-Z]' OUTPUT_VARIABLE dirs)
    string(LENGTH ${dirs} string_length)
    string(SUBSTRING ${dirs} 0 ${string_length}-2 dirs_trimmed)
    string(REPLACE "\n" ";" dirs_list ${dirs_trimmed})

    set(lib_dirs_list)
    set(lib_prod_list)
    foreach(prod ${dirs_list})
#message(STATUS "test ${prod}")
#message(STATUS "var $ENV{${prod}_LIB}")
        if(prod STREQUAL "")
            continue()
        endif(prod STREQUAL "")
        if(IS_DIRECTORY $ENV{${prod}_LIB})
#message(STATUS "$ENV{${prod}_LIB}")
            list(APPEND lib_prod_list ${prod})
            list(APPEND lib_dirs_list $ENV{${prod}_LIB})
            insert_into_map(${prod} $ENV{${prod}_LIB})
        endif(IS_DIRECTORY $ENV{${prod}_LIB})
    endforeach(prod)
#message(STATUS "${lib_dirs_list}")

endmacro()
