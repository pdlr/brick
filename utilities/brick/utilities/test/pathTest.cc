/**
***************************************************************************
* @file brick/utilities/test/pathTest.cc
*
* Source file defining tests for file handling.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <string>

#include <brick/utilities/path.hh>
#include <brick/test/testFixture.hh>


// We disable a bunch of tests due to unsafe file IO.
#ifndef BRICK_UNITTEST_USE_UNSAFE_FILEIO
#define BRICK_UNITTEST_USE_UNSAFE_FILEIO 0
#endif

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class PathTest : public TestFixture<PathTest> {

    public:

      PathTest();
      ~PathTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testIsDirectory();
      void testJoinPath();
      void testListDirectory();
      void testRecursiveListDirectory();
      void testSplitExtension();
      void testSplitPath();

    }; // class PathTest


    /* ============== Member Function Definititions ============== */

    PathTest::
    PathTest()
      : TestFixture<PathTest>("PathTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testIsDirectory);
      BRICK_TEST_REGISTER_MEMBER(testJoinPath);
      BRICK_TEST_REGISTER_MEMBER(testListDirectory);
      BRICK_TEST_REGISTER_MEMBER(testRecursiveListDirectory);
      BRICK_TEST_REGISTER_MEMBER(testSplitExtension);
      BRICK_TEST_REGISTER_MEMBER(testSplitPath);
    }


    void
    PathTest::
    testIsDirectory()
    {

#if BRICK_UNITTEST_USE_UNSAFE_FILEIO

      // Set up a test directory.
      system("rm -rf _disposable_test_directory_");
      system("mkdir _disposable_test_directory_");
      system("touch _disposable_test_directory_/foo.txt");
      system("touch _disposable_test_directory_/bar.png");
      system("touch _disposable_test_directory_/baz.sxw");

      // Verify results.
      try{
        BRICK_TEST_ASSERT(isDirectory("_disposable_test_directory_"));
        BRICK_TEST_ASSERT(!isDirectory("_disposable_test_directory_/foo.txt"));
        BRICK_TEST_ASSERT(!isDirectory("_disposable_test_directory_/bar.png"));
        BRICK_TEST_ASSERT(!isDirectory("_disposable_test_directory_/baz.sxw"));
      } catch(const TestException&) {
        // Clean up.
        system("rm -rf _disposable_test_directory_");
        throw;
      }

      // Clean up.
      system("rm -rf _disposable_test_directory_");

#endif /* #if BRICK_UNITTEST_USE_UNSAFE_FILEIO */

    }


    void
    PathTest::
    testJoinPath()
    {
#ifdef _WIN32
      BRICK_TEST_ASSERT(joinPath("foo", "bar.baz") == "foo\\bar.baz");
      BRICK_TEST_ASSERT(joinPath("foo\\bar", "baz.hum") == "foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo\\bar", "") == "foo\\bar");
      BRICK_TEST_ASSERT(joinPath("\\foo", "bar\\baz.hum") == "\\foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("", "bar\\baz.hum") == "bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("", "\\bar\\baz.hum") == "\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\foo", "bar\\baz.hum") == "\\foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\foo", "\\bar\\baz.hum") == "\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo\\bar", "\\baz.hum") == "foo\\bar\\\\baz.hum");

      BRICK_TEST_ASSERT(joinPath("foo\\", "bar.baz") == "foo\\bar.baz");
      BRICK_TEST_ASSERT(joinPath("foo\\bar\\", "baz.hum") == "foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo\\bar\\", "") == "foo\\bar\\");
      BRICK_TEST_ASSERT(joinPath("\\foo\\", "bar\\baz.hum")
                        == "\\foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\", "bar\\baz.hum") == "\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\", "\\bar\\baz.hum") == "\\\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\foo\\", "bar\\baz.hum")
                        == "\\foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("\\foo\\", "\\bar\\baz.hum")
                        == "\\bar\\baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo\\bar\\", "\\baz.hum")
                        == "foo\\bar\\\\baz.hum");

#if 0
      // Hmm.  Not sure what the right behavior is here.
      // Should joinPath("foo", "/bar") give "foo//bar", "foo/bar", or "/bar"?

      // Test cases suggested by users.
      BRICK_TEST_ASSERT(joinPath("testVariety", "\\training\\phase1model\\")
                        == "testVariety\\training\\phase1model\\");
      BRICK_TEST_ASSERT(joinPath("\\testVariety", "\\training\\phase1model\\")
                        == "\\testVariety\\training\\phase1model\\");
      BRICK_TEST_ASSERT(joinPath("testVariety\\", "\\training\\phase1model\\")
                        == "testVariety\\training\\phase1model\\");
      BRICK_TEST_ASSERT(joinPath("\\testVariety\\", "\\training\\phase1model\\")
                        == "\\testVariety\\training\\phase1model\\");
#endif /* #if 0 */

#else
      BRICK_TEST_ASSERT(joinPath("foo", "bar.baz") == "foo/bar.baz");
      BRICK_TEST_ASSERT(joinPath("foo/bar", "baz.hum") == "foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo/bar", "") == "foo/bar");
      BRICK_TEST_ASSERT(joinPath("/foo", "bar/baz.hum") == "/foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("", "bar/baz.hum") == "bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("", "/bar/baz.hum") == "/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/foo", "bar/baz.hum") == "/foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/foo", "/bar/baz.hum") == "/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo/bar", "/baz.hum") == "/baz.hum");

      BRICK_TEST_ASSERT(joinPath("foo/", "bar.baz") == "foo/bar.baz");
      BRICK_TEST_ASSERT(joinPath("foo/bar/", "baz.hum") == "foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo/bar/", "") == "foo/bar/");
      BRICK_TEST_ASSERT(joinPath("/foo/", "bar/baz.hum") == "/foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/", "bar/baz.hum") == "/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/", "/bar/baz.hum") == "/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/foo/", "bar/baz.hum") == "/foo/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("/foo/", "/bar/baz.hum") == "/bar/baz.hum");
      BRICK_TEST_ASSERT(joinPath("foo/bar/", "/baz.hum") == "/baz.hum");

#if 0
      // Hmm.  Not sure what the right behavior is here.
      // Should joinPath("foo", "/bar") give "foo//bar", "foo/bar", or "/bar"?

      // Test cases suggested by users.
      BRICK_TEST_ASSERT(joinPath("testVariety", "/training/phase1model/")
                        == "/training/phase1model/");
      BRICK_TEST_ASSERT(joinPath("/testVariety", "/training/phase1model/")
                        == "/training/phase1model/");
      BRICK_TEST_ASSERT(joinPath("testVariety/", "/training/phase1model/")
                        == "/training/phase1model/");
      BRICK_TEST_ASSERT(joinPath("/testVariety/", "/training/phase1model/")
                        == "/training/phase1model/");
      BRICK_TEST_ASSERT(false);   // One of these new tests should fail,
      // or else talk with Stager about his
      // email of 10/1/09.  Remove this line
      // once the bug is fixed.
#endif /* #if 0 */

#endif /* #ifdef _WIN32 */
    }


    void
    PathTest::
    testListDirectory()
    {

#if BRICK_UNITTEST_USE_UNSAFE_FILEIO

      // Set up a test directory.
      system("rm -rf _disposable_test_directory_");
      system("mkdir _disposable_test_directory_");
      system("touch _disposable_test_directory_/foo.txt");
      system("touch _disposable_test_directory_/bar.png");
      system("touch _disposable_test_directory_/baz.sxw");

      // List the contents and sort them into a predictable order.
      std::vector<std::string> listing =
        listDirectory("_disposable_test_directory_");
      std::sort(listing.begin(), listing.end());

      // Verify results.
      try{
        BRICK_TEST_ASSERT(listing.size() == 5);
        BRICK_TEST_ASSERT(listing[0] == ".");
        BRICK_TEST_ASSERT(listing[1] == "..");
        BRICK_TEST_ASSERT(listing[2] == "bar.png");
        BRICK_TEST_ASSERT(listing[3] == "baz.sxw");
        BRICK_TEST_ASSERT(listing[4] == "foo.txt");
      } catch(const TestException&) {
        // Clean up.
        system("rm -rf _disposable_test_directory_");
        throw;
      }

      // List the contents again and sort them into a predictable order.
      // This time ask for full path.
      listing =
        listDirectory("_disposable_test_directory_", true);
      std::sort(listing.begin(), listing.end());

      // Verify results.
      try{
        BRICK_TEST_ASSERT(listing.size() == 5);
        BRICK_TEST_ASSERT(listing[0] == "_disposable_test_directory_/.");
        BRICK_TEST_ASSERT(listing[1] == "_disposable_test_directory_/..");
        BRICK_TEST_ASSERT(listing[2] == "_disposable_test_directory_/bar.png");
        BRICK_TEST_ASSERT(listing[3] == "_disposable_test_directory_/baz.sxw");
        BRICK_TEST_ASSERT(listing[4] == "_disposable_test_directory_/foo.txt");
      } catch(const TestException&) {
        // Clean up.
        system("rm -rf _disposable_test_directory_");
        throw;
      }

      // Clean up.
      system("rm -rf _disposable_test_directory_");

#endif /* #if BRICK_UNITTEST_USE_UNSAFE_FILEIO */

    }


    void
    PathTest::
    testRecursiveListDirectory()
    {

#if BRICK_UNITTEST_USE_UNSAFE_FILEIO

      // Set up a test directory.
      system("rm -rf _disposable_test_directory_");
      system("mkdir _disposable_test_directory_");
      system("touch _disposable_test_directory_/foo.txt");
      system("touch _disposable_test_directory_/hoo.sxw");
      system("mkdir _disposable_test_directory_/subdirectory0");
      system("touch _disposable_test_directory_/subdirectory0/bar.png");
      system("touch _disposable_test_directory_/subdirectory0/baz.png");
      system("mkdir _disposable_test_directory_/subdirectory0/subdirectory1");

      // List the contents and sort them into a predictable order.
      std::vector<std::string> listing =
        recursiveListDirectory("_disposable_test_directory_");
      std::sort(listing.begin(), listing.end());

      // Verify results.
      try{
        size_t index0 = 0;
        BRICK_TEST_ASSERT(listing.size() == 4);
        BRICK_TEST_ASSERT(listing[index0++] == "foo.txt");
        BRICK_TEST_ASSERT(listing[index0++] == "hoo.sxw");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/bar.png");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/baz.png");
      } catch(const TestException&) {
        // Clean up.
        // system("rm -rf _disposable_test_directory_");
        throw;
      }

      // List the contents and sort them into a predictable order.  This
      // time ask for directory names to be included.
      listing =
        recursiveListDirectory("_disposable_test_directory_", false, true);
      std::sort(listing.begin(), listing.end());

      // Verify results.
      try{
        size_t index0 = 0;
        BRICK_TEST_ASSERT(listing.size() == 12);
        BRICK_TEST_ASSERT(listing[index0++] == ".");
        BRICK_TEST_ASSERT(listing[index0++] == "..");
        BRICK_TEST_ASSERT(listing[index0++] == "foo.txt");
        BRICK_TEST_ASSERT(listing[index0++] == "hoo.sxw");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/.");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/..");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/bar.png");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/baz.png");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/subdirectory1");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/subdirectory1/.");
        BRICK_TEST_ASSERT(listing[index0++] == "subdirectory0/subdirectory1/..");
      } catch(const TestException&) {
        // Clean up.
        system("rm -rf _disposable_test_directory_");
        throw;
      }

      // List the contents and sort them into a predictable order.  This
      // time ask for directory names to be included and request full
      // path name.
      listing =
        recursiveListDirectory("_disposable_test_directory_", true, true);
      std::sort(listing.begin(), listing.end());

      // Verify results.
      try{
        size_t index0 = 0;
        BRICK_TEST_ASSERT(listing.size() == 12);
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/.");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/..");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/foo.txt");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/hoo.sxw");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/subdirectory0");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/subdirectory0/.");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/subdirectory0/..");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/subdirectory0/bar.png");
        BRICK_TEST_ASSERT(listing[index0++]
                          == "_disposable_test_directory_/subdirectory0/baz.png");
        BRICK_TEST_ASSERT(
          listing[index0++]
          == "_disposable_test_directory_/subdirectory0/subdirectory1");
        BRICK_TEST_ASSERT(
          listing[index0++]
          == "_disposable_test_directory_/subdirectory0/subdirectory1/.");
        BRICK_TEST_ASSERT(
          listing[index0++]
          == "_disposable_test_directory_/subdirectory0/subdirectory1/..");
      } catch(const TestException&) {
        // Clean up.
        system("rm -rf _disposable_test_directory_");
        throw;
      }

      // Clean up.
      system("rm -rf _disposable_test_directory_");

#endif /* #if BRICK_UNITTEST_USE_UNSAFE_FILEIO */

    }


    void
    PathTest::
    testSplitExtension()
    {
#ifdef _WIN32
      // Try splitting a normal looking filename.
      std::pair<std::string, std::string> basename_extension =
        splitExtension("\\foo\\bar.baz");
      BRICK_TEST_ASSERT(basename_extension.first == "\\foo\\bar");
      BRICK_TEST_ASSERT(basename_extension.second == ".baz");

      // Make sure we don't get faked out by non-extension '.'
      // characters.
      basename_extension = splitExtension("\\foo\\bar.j\\baz");
      BRICK_TEST_ASSERT(basename_extension.first == "\\foo\\bar.j\\baz");
      BRICK_TEST_ASSERT(basename_extension.second == "");

      // Check with a filename that's nothing but an extension.
      basename_extension = splitExtension(".baz");
      BRICK_TEST_ASSERT(basename_extension.first == "");
      BRICK_TEST_ASSERT(basename_extension.second == ".baz");
#else
      // Try splitting a normal looking filename.
      std::pair<std::string, std::string> basename_extension =
        splitExtension("/foo/bar.baz");
      BRICK_TEST_ASSERT(basename_extension.first == "/foo/bar");
      BRICK_TEST_ASSERT(basename_extension.second == ".baz");

      // Make sure we don't get faked out by non-extension '.'
      // characters.
      basename_extension = splitExtension("/foo/bar.j/baz");
      BRICK_TEST_ASSERT(basename_extension.first == "/foo/bar.j/baz");
      BRICK_TEST_ASSERT(basename_extension.second == "");

      // Check with a filename that's nothing but an extension.
      basename_extension = splitExtension(".baz");
      BRICK_TEST_ASSERT(basename_extension.first == "");
      BRICK_TEST_ASSERT(basename_extension.second == ".baz");
#endif
    }

    void
    PathTest::
    testSplitPath()
    {
      std::pair<std::string, std::string> directoryName_fileName;

#ifdef _WIN32
      directoryName_fileName = splitPath("\\foo\\bar\\baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "\\foo\\bar\\");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("bar\\baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "bar\\");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("\\foo\\bar\\baz.hum\\");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "\\foo\\bar\\baz.hum\\");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "");
#else
      directoryName_fileName = splitPath("/foo/bar/baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "/foo/bar/");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("bar/baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "bar/");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("baz.hum");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "baz.hum");

      directoryName_fileName = splitPath("/foo/bar/baz.hum/");
      BRICK_TEST_ASSERT(directoryName_fileName.first == "/foo/bar/baz.hum/");
      BRICK_TEST_ASSERT(directoryName_fileName.second == "");
#endif
    }

  } // namespace utilities

} // namespace brick


#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::PathTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::PathTest currentTest;

}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
