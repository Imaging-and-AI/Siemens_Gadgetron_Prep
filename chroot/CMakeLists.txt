set(ENV_NAME "$ENV{ENV_NAME}")

configure_file("start-conda-gadgetron.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-conda-gadgetron.sh @ONLY)

install(PROGRAMS 
  ${CMAKE_CURRENT_BINARY_DIR}/start-conda-gadgetron.sh
  DESTINATION ${GADGETRON_INSTALL_CHROOT_SCRIPTS_PATH}
  COMPONENT scripts
)
