# Common Ambient Variables:
#   CURRENT_BUILDTREES_DIR    = ${VCPKG_ROOT_DIR}\buildtrees\${PORT}
#   CURRENT_PACKAGES_DIR      = ${VCPKG_ROOT_DIR}\packages\${PORT}_${TARGET_TRIPLET}
#   CURRENT_PORT_DIR          = ${VCPKG_ROOT_DIR}\ports\${PORT}
#   PORT                      = current port name (zlib, etc)
#   TARGET_TRIPLET            = current triplet (x86-windows, x64-windows-static, etc)
#   VCPKG_CRT_LINKAGE         = C runtime linkage type (static, dynamic)
#   VCPKG_LIBRARY_LINKAGE     = target library linkage type (static, dynamic)
#   VCPKG_ROOT_DIR            = <C:\path\to\current\vcpkg>
#   VCPKG_TARGET_ARCHITECTURE = target architecture (x64, x86, arm)
#

include(vcpkg_common_functions)

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO mfem/mfem
	SHA512 e79a0787ac125d887fdb02b69ea58c6182a32a300718a2ff9b09adb2993a0071c7d2aa31f2b0553e6305750b99f6e041190e8945bd20b88eb5c5ae46b8101de9
    REF master
)

#set(SOURCE_PATH ${CURRENT_BUILDTREES_DIR}/src)

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA # Disable this option if project cannot be built with Ninja
	OPTIONS -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=1
    # OPTIONS -DUSE_THIS_IN_ALL_BUILDS=1 -DUSE_THIS_TOO=2
    #OPTIONS_RELEASE -DOPTIMIZE=1
    #OPTIONS_DEBUG -DDEBUGGABLE=1
)

vcpkg_install_cmake()

file(COPY ${CURRENT_PACKAGES_DIR}/lib/mfem.dll DESTINATION ${CURRENT_PACKAGES_DIR}/bin)
file(COPY ${CURRENT_PACKAGES_DIR}/debug/lib/mfem.dll DESTINATION ${CURRENT_PACKAGES_DIR}/debug/bin)  
file(REMOVE ${CURRENT_PACKAGES_DIR}/lib/mfem.dll ${CURRENT_PACKAGES_DIR}/debug/lib/mfem.dll)
file(COPY ${CURRENT_PACKAGES_DIR}/lib/cmake/mfem DESTINATION ${CURRENT_PACKAGES_DIR}/share)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/lib/cmake)
file(COPY ${CURRENT_PACKAGES_DIR}/debug/lib/cmake/mfem/MFEMTargets-debug.cmake DESTINATION ${CURRENT_PACKAGES_DIR}/share)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/lib/cmake)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

# Handle copyright
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/mfem RENAME copyright)

vcpkg_copy_pdbs()

file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)

