function(import_delib)
  if(NOT TARGET delib::Delib)
    include(DownloadProject)
    download_project(
      PROJ delib
      GIT_REPOSITORY      git@gitlab.chab.ethz.ch:reiher/delib.git
      GIT_TAG             v0.2.1
      UPDATE_DISCONNECTED 1
    )
    add_subdirectory(${delib_SOURCE_DIR} ${delib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endfunction()
