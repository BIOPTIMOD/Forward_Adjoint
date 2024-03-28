
add_library(adj STATIC
   bioptimod_memory.f90 fm34.f90 SLAE.f90 adj_new.f90
)

add_executable(adj_exe
   main_adj.f90
)
set_property(TARGET adj_exe PROPERTY OUTPUT_NAME adj)

target_link_libraries(adj_exe PRIVATE adj)

