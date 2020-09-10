function(download_file source_file target_file)
  if(NOT EXISTS ${target_file})
    file(DOWNLOAD ${source_file} ${target_file} STATUS PULL_STATUS)
    list(GET PULL_STATUS 0 PULL_RESULT)
    if(NOT ${PULL_RESULT} EQUAL 0)
      if(EXISTS ${target_file})
        file(REMOVE ${target_file})
      endif()
      list(GET PULL_STATUS 1 PULL_ERROR)
      cmessage(FATAL_ERROR
        "Could not download file: ${PULL_ERROR}"
        "Please download the file at ${source_file}"
        "and place it in your build tree at ${target_file}."
      )
    endif()
  endif()
endfunction()
