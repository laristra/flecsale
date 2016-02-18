# Summary

Individual user apps that test the functionality of the ALE library go
here.  Apps are essentially the main drivers that test and execute
different functionallity developped in `<PROJECT_ROOT>/src`

To add an app, create a translation unit ("*.cc" file) with a main
function in this directory and add it to the `CMakeLists.txt` file:

```cmake
add_executable( dummy_app dummy_app.cc )
```

If you need to link against any libraries, the main **ale** project
library for example, you need to add the following line next:

```cmake
target_link_libraries( dummy_app ale )
```

Building the project will create an executable `dummy_app` in
`build\apps`.  This **dummy** example is an actual example created
within this project to help you understand.
