set(SRC
    src/uni_bis.c
)

add_library(opm ${SRC})
target_include_directories(opm PUBLIC include/)
target_link_libraries(opm PUBLIC stdx)