# Summary

This is where the main ALE library functionallity is developped.  The
end result of the code whithin this folder is a library,
**libale.a**, that users can link against in order to create actual
applications.  There is no `main` function defined anywhere within
this directory.

Developpers should create and test their functionality here.  A sample
sub-package is provided for instructional purposes in `dummy`.
Developpers should familiarize themselves with the
[Cinch](https://github.com/losalamos/cinch) build system and Google's
C++ testing framework
[gtest](https://github.com/google/googletest/blob/master/googletest/docs/Primer.md).

Some key points to remember, all code within `src` is protected within
the `ale` namespace, and each sub-folder/sub-package has its own
namespace.  The latter helps help group functionality together and
make it easy to tell where functions are defined.  Additionally tests
and test-specific code goes inside folders named `test` that live where
the main functionality was defined.

