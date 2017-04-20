/**
***************************************************************************
* @file brick/utilities/path.hh
*
* Header file declaring routines for working with the filesystem
*
* Copyright (C) 2003-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_PATH_HH
#define BRICK_UTILITIES_PATH_HH

#include <string>
#include <vector>

namespace brick {

  namespace utilities {

    /** 
     * This function returns a bool indicating whether or not the
     * specified path is a directory.  If the path is a symbolic link,
     * the return value is currently unspecified, but will eventually be
     * true iff the link points to a directory.
     * 
     * @param path This argument is the filename to evaluate.
     * 
     * @return The return value is true if path refers to a directory,
     * false otherwise.
     */
    bool
    isDirectory(std::string const& path);
  
  
    /** 
     * This function returns a bool indicating whether or not the
     * specified file or directory exists.  If the path is a symbolic link,
     * the return value is currently unspecified, but will eventually be
     * true.
     * 
     * @param path This argument is the filename to evaluate.
     * 
     * @return The return value is true if path exists, false
     * otherwise.
     */
    bool
    isExistingPath(std::string const& path);
  

    /** 
     * This function reolves references to ".", "..", symbolic links,
     * etc., and returns the canonicalized absolute pathname
     * corresponding to its input.
     * 
     * @param inputPath The path to be canonicalized.
     * 
     * @return The return value is the resolved path, or an empty
     * string if path resolution fails.
     */
    std::string
    getAbsolutePath(std::string const& inputPath);
    

    /** 
     * This function returns a bool indicating whether or not the
     * specified file exists and is not a directory or other special
     * file.  If the path is a symbolic link, the return value is
     * currently unspecified, but will eventually be true iff the link
     * points to a file.
     * 
     * @param path This argument is the filename to evaluate.
     * 
     * @return The return value is true if path refers to a regular
     * file, false otherwise.
     */
    bool
    isRegularFile(std::string const& path);
    
    
    /**
     * Joins two path elements with the appropriate path delimiter.
     * For example:
     *
     *   joinPath("foo", "bar");
     * 
     * might give
     *
     *   "foo/bar"
     * 
     * while
     * 
     *   joinPath("foo/baz/", "bar");
     * 
     * might give
     * 
     *   "/foo/baz/bar"
     **/
    std::string
    joinPath(std::string const& part0, std::string const& part1);

  
    /**
     * Returns the names of the entries in the specified directory, in
     * no particular order.  For example,
     * 
     *   listDirectory("/etc");
     * 
     * might give
     * 
     *   ["fstab", "init.d", "modules.conf", ...]
     *
     * while
     *
     *   listDirectory("/etc", true);
     * 
     * might give
     * 
     *   ["/etc/fstab", "/etc/init.d", "/etc/modules.conf", ...]
     **/
    std::vector<std::string>
    listDirectory(std::string const& directoryName, bool fullPath=false);


    /** 
     * Returns the names of files in the directory tree below the
     * specified directory.  For example,
     *
     *   recursiveListDirectory("/etc");
     *
     * might give
     *
     *   ["fstab", "init.d/cron", "init.d/ssh", "modules.conf", ...]
     *
     * while
     *
     *   recursiveListDirectory("/etc", false, true);
     *
     * might give
     *
     *   ["fstab", "init.d", "init.d/cron", "init.d/ssh", "modules.conf", ...]
     *
     * and
     *
     *   recursiveListDirectory("/etc", true, true);
     *
     * might give
     *
     *   ["/etc/fstab", "/etc/init.d", "/etc/init.d/cron", "/etc/init.d/ssh",
     *    "/etc/modules.conf", ...]
     *
     * @param directoryName This argument specifies the directory to be listed.
     *
     * @param fullPath Indicates whether the returned filenames should
     * be prepended with directoryName.
     *
     * @param includeDirectoryNames Indicates whether the names of
     * subdirectories should be included in the listing.
     *
     * @return The return value is a vector of file names.
     */
    std::vector<std::string>
    recursiveListDirectory(std::string const& directoryName,
                           bool fullPath=false,
                           bool includeDirectoryNames=false);


    /** 
     * This function searches a sequence of directories looking for
     * the specified file.  If the file is found, the complete path to
     * the found file is returned through reference argument fullPath.
     * If argument discardInputDir is true, then any preceding
     * directory information will be stripped from the input fileName
     * before searching.  If the requested file exists in more than
     * one directory, the first match will be returned.
     *
     * For example, imagine the following directory structure:
     *
     * @code
     *   /home
     *     /rick
     *       /dir0
     *         /subdir0
     *           foo.txt
     *       /dir1
     *         bar.txt
     *       /dir2
     *         foo.txt
     *       /dir3
     *         bar.txt
     * @endcode
     *
     * And assume the following vector of directory names:
     *
     * @code
     *   std::vector<std::string> searchPath(3);
     *   searchPath[0] = /home/rick/dir0;
     *   searchPath[1] = /home/rick/dir1;
     *   searchPath[2] = /home/rick/dir2;
     * @endcode
     *
     * With this directory structure, the call:
     *
     * @code
     *   std::string result;
     *   searchForFile("subdir0/foo.txt", searchPath.begin(), searchPath.end(),
     *                 result, false);
     * @endcode
     *
     * would return true and set result to
     * "/home/rick/dir0/subdir0/foo.txt".
     *
     * The call:
     * 
     * @code
     *   std::string result;
     *   searchForFile("subdir0/foo.txt", searchPath.begin(), searchPath.end(),
     *                 result, true);
     * @endcode
     *
     * would return true and set result to
     * "/home/rick/dir2/foo.txt".
     *
     * The call:
     *
     * @code
     *   std::string result;
     *   searchForFile("subdir0/baz.txt", searchPath.begin(), searchPath.end(),
     *                 result, true);
     * @endcode
     *
     * would return false and not touch result.
     *   
     * @param fileName This argument is the file to be found.
     * 
     * @param pathBegin This argument is an interator pointing to the
     * beginning of a sequence of directories in which to search.
     * 
     * @param pathEnd This argument is an interator pointing to the
     * end (in the STL sense) of a sequence of directories in which to
     * search.
     * 
     * @param fullPath This argument is used to return the full path
     * to the located file.  If the file is not found, this argument
     * will not be accessed.
     *
     * @param discardInputDir If this argument is set to true, than
     * any preceding directory information will be removed from
     * argument fileName before the search is conducted.  See the
     * examples above for more information.
     * 
     * @return The return value is true if the file is found, false otherwise.
     */
    template <class IterType>
    bool
    searchForFile(std::string const& fileName,
                  IterType pathBegin, IterType pathEnd,
                  std::string& fullPath,
                  bool discardInputDir = false);
    
  
    /**
     * Returns a std::pair<std::string, std::string> containing the fileName
     * without its extension, and the extension.  For example:
     *
     *   splitExt("/foo/bar.baz")
     * 
     * returns
     *
     *   {"/foo/bar", ".baz"}
     *
     * @param filename This is the filename to be split.
     **/
    std::pair<std::string, std::string>
    splitExtension(std::string const& fileName);


    /** 
     * This function accepts a path returns a pair of strings in which
     * the first element is the directory name and the second is the
     * filename.  For example:
     *
     *   splitPath("/foo/bar.baz")
     *
     * might return
     *
     *   {"/foo/", "bar.baz"}
     *
     * Also,
     *
     *   splitPath("bar.baz")
     *
     * might return
     *
     *   {"", "bar.baz"}
     *
     *   splitPath("/foo/")
     *
     * might return
     *
     *   {"/foo/", ""}
     *
     * @param path This is the pathname to be split.
     * 
     * @return The return value is a pair containing first the directory
     * name, and second the file name.
     */
    std::pair<std::string, std::string>
    splitPath(std::string const& path);

  } // namespace utilities
    
} // namespace brick


/* ======= Declarations of inline and template functions. ======= */

namespace brick {

  namespace utilities {
    
    // This function searches a sequence of directories looking for
    // the specified file.
    template <class IterType>
    bool
    searchForFile(std::string const& fileName,
                  IterType pathBegin, IterType pathEnd,
                  std::string& fullPath,
                  bool discardInputDir)
    {
      std::string baseName;
      if(discardInputDir) {
        baseName = splitPath(fileName).second;
      } else {
        baseName = fileName;
      }
      while(pathBegin != pathEnd) {
        std::string candidateName = joinPath(*pathBegin, baseName);
        if(isRegularFile(candidateName)) {
          fullPath = candidateName;
          return true;
        }
        ++pathBegin;
      }
      return false;
    }

  } // namespace utilities
  
} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_PATH_HH */
