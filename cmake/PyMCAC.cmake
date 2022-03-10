find_package(PythonInterp 3.8 REQUIRED)

# Generate the venv
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${PYTHON_EXECUTABLE} -m venv ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -U pip setuptools wheel
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -U -r
            ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
    COMMENT "Creating a venv with PyMCAC dependencies"
)

# install/update pymcac in the venv
add_custom_target(
    pymcac ALL
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/venv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/pip install -e ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Installing PyMCAC in the venv"
)
