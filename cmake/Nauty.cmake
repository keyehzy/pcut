# Private, pinned implementation dependency. Configure upstream headers in the
# binary tree; the downloaded source tree is never patched or configured in place.
include(FetchContent)
FetchContent_Declare(nauty
  URL https://pallini.di.uniroma1.it/nauty2_9_3.tar.gz
  URL_HASH SHA256=9fc4edae04f88a0f5883985be3b39cf7f898fd6cc96e96b9ee25452743cc1b5b)
FetchContent_MakeAvailable(nauty)
file(MAKE_DIRECTORY "${nauty_BINARY_DIR}")
string(SHA256 nauty_configuration "${nauty_SOURCE_DIR};${CMAKE_C_COMPILER};${CMAKE_C_FLAGS};tls;wordsize64")
set(nauty_configured "")
if(EXISTS "${nauty_BINARY_DIR}/pcut-configuration")
  file(READ "${nauty_BINARY_DIR}/pcut-configuration" nauty_configured)
endif()
if(NOT nauty_configured STREQUAL nauty_configuration OR NOT EXISTS "${nauty_BINARY_DIR}/nauty.h")
  execute_process(COMMAND sh "${nauty_SOURCE_DIR}/configure"
    "--srcdir=${nauty_SOURCE_DIR}" --enable-tls --enable-wordsize=64
    "CC=${CMAKE_C_COMPILER}" "CFLAGS=${CMAKE_C_FLAGS}"
    WORKING_DIRECTORY "${nauty_BINARY_DIR}"
    OUTPUT_FILE "${nauty_BINARY_DIR}/configure.log"
    ERROR_FILE "${nauty_BINARY_DIR}/configure-errors.log"
    RESULT_VARIABLE nauty_configure_result)
  if(NOT nauty_configure_result EQUAL 0)
    message(FATAL_ERROR "nauty configure failed: ${nauty_BINARY_DIR}/configure-errors.log")
  endif()
  file(WRITE "${nauty_BINARY_DIR}/pcut-configuration" "${nauty_configuration}")
endif()
add_library(pcut_nauty STATIC
  ${nauty_SOURCE_DIR}/nauty.c ${nauty_SOURCE_DIR}/nautil.c
  ${nauty_SOURCE_DIR}/nausparse.c ${nauty_SOURCE_DIR}/naugraph.c
  ${nauty_SOURCE_DIR}/schreier.c ${nauty_SOURCE_DIR}/naurng.c
  ${PROJECT_SOURCE_DIR}/src/nauty_alloc.cpp)
set_target_properties(pcut_nauty PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(pcut_nauty PRIVATE c_std_11 cxx_std_20)
target_include_directories(pcut_nauty PUBLIC
  $<BUILD_INTERFACE:${nauty_BINARY_DIR}> $<BUILD_INTERFACE:${nauty_SOURCE_DIR}>)
if(NOT CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
  message(FATAL_ERROR "The nauty allocation adapter currently requires GCC or Clang")
endif()
target_compile_options(pcut_nauty PRIVATE
  $<$<COMPILE_LANGUAGE:C>:-fexceptions;-include;${PROJECT_SOURCE_DIR}/src/nauty_alloc.h>)
if(PCUT_SANITIZE)
  target_compile_options(pcut_nauty PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(pcut_nauty PUBLIC -fsanitize=address,undefined)
endif()
install(TARGETS pcut_nauty EXPORT pcutTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${nauty_SOURCE_DIR}/COPYRIGHT ${nauty_SOURCE_DIR}/LICENSE-2.0.txt
  DESTINATION ${CMAKE_INSTALL_DATADIR}/pcut/nauty)
