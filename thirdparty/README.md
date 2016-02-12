# Obtaining the Thirdparty Libraries
If you already have the necessary libraries in a folder somewhere, you can skip this step.  Otherwise, the build system can download them for you.  A good place to put them is inside the reposotory itself.  An empty folder was created inside the git repository for this purpose.  It is located in **<ALE_DIR>/thirdparty/files** where **<ALE_DIR>** is the location of the cloned git repository.

    TPL_DOWNLOAD_PATH=<ALE_DIR>/thirdparty/files <ALE_DIR>/arch/download-tpl.sh

