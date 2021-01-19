

find_package(PythonInterp 3.6 REQUIRED)

# Generate the virtualenv
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${PYTHON_EXECUTABLE} -m venv ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -U pip setuptools wheel
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -U numpy
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
)

# install/update pymcac in the venv
add_custom_target(pymcac
    ALL
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -e ${CMAKE_CURRENT_SOURCE_DIR}
)
