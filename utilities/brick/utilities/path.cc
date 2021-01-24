/**
***************************************************************************
* @file brick/utilities/path.cc
*
* Source file defining routines for working with the filesystem
*
* Copyright (C) 2003-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include <list>
#include <brick/portability/filesystem.hh>
#include <brick/utilities/path.hh>

namespace brick {

  namespace utilities {

    using brick::portability::pathDelimiter;
    using brick::portability::extensionDelimiter;


    // This function returns a bool indicating whether or not the
    // specified path is a directory.
    bool
    isDirectory(std::string const& path)
    {
      return brick::portability::isDirectory(path);
    }


    bool
    isExistingPath(std::string const& path)
    {
      struct stat statBuf;
      memset(&statBuf, 0, sizeof(statBuf));
      int returnValue = stat(path.c_str(), &statBuf);
      if(returnValue == 0) {
        return true;
      }
      return false;
    }


    // This function returns a bool indicating whether or not the
    // specified path is a regular file.
    bool
    isRegularFile(std::string const& path)
    {
      struct stat statBuf;
      memset(&statBuf, 0, sizeof(statBuf));
      int returnValue = stat(path.c_str(), &statBuf);
      if(returnValue == 0) {
        // if(S_ISREG(statBuf.st_mode)) {
        if((statBuf.st_mode & S_IFREG) == S_IFREG) {
          return true;
        }
      }
      return false;
    }


    // This function reolves references to ".", "..", symbolic links,
    // etc., and returns the canonicalized absolute pathname
    // corresponding to its input.
    std::string
    getAbsolutePath(std::string const& inputPath)
    {
      return brick::portability::getAbsolutePath(inputPath);
    }


    // Joins two path elements with the appropriate delimiter.
    std::string
    joinPath(std::string const& part0, std::string const& part1)
    {
      return brick::portability::joinPath(part0, part1);
    }


    // Returns the names of the entries in the specified directory, in
    // no particular order.
    std::vector<std::string>
    listDirectory(std::string const& directoryName, bool fullPath)
    {
      // Dispatch to dlrPortability function.
      return brick::portability::listDirectory(directoryName, fullPath);
    }


    // Returns the names of files in the directory tree below the
    // specified directory.
    std::vector<std::string>
    recursiveListDirectory(std::string const& directoryName,
                           bool fullPath,
                           bool includeDirectoryNames)
    {
      typedef std::vector<std::string>::const_iterator FileNameIter;

      // We'll accumulate the total listing in fileNameList.
      std::list<std::string> fileNameList;

      // For starters, just find the entries in the specified directory.
      std::vector<std::string> topLevelFileNameList =
        listDirectory(directoryName);

      // Process each entry in turn.
      for(FileNameIter fileNameIter = topLevelFileNameList.begin();
          fileNameIter != topLevelFileNameList.end();
          ++fileNameIter) {

        std::string fullName = joinPath(directoryName, *fileNameIter);
        if(isDirectory(fullName)) {
          // Looks like this entry is a directory.  Add it, if
          // appropriate.
          if(includeDirectoryNames) {
            fileNameList.push_back(*fileNameIter);
          }

          if(*fileNameIter != "." && *fileNameIter != "..") {
            // Recurse into non-trivial directory entries, and add the
            // result to our accumulated list.
            std::vector<std::string> subListing =
              recursiveListDirectory(fullName, false, includeDirectoryNames);
            std::vector<std::string>::const_iterator iter = subListing.begin();
            while(iter != subListing.end()) {
              fileNameList.push_back(joinPath(*fileNameIter, *iter));
              ++iter;
            }
          }
        } else {  // if(isDirectory(...))
          fileNameList.push_back(*fileNameIter);
        }
      }

      // Now copy the file names into a vector.  Adding the full path,
      // if appriate.
      std::vector<std::string> finalVector(fileNameList.size());
      if(fullPath) {
        std::list<std::string>::const_iterator iter = fileNameList.begin();
        std::vector<std::string>::iterator finalIter = finalVector.begin();
        while(iter != fileNameList.end()) {
          *finalIter = joinPath(directoryName, *iter);
          ++finalIter;
          ++iter;
        }
      } else {
        std::copy(fileNameList.begin(), fileNameList.end(),
                  finalVector.begin());
      }
      return finalVector;
    }


    // Returns a std::pair<std::string, std::string> containing the fileName
    // without its extension, and the extension.
    std::pair<std::string, std::string>
    splitExtension(std::string const& fileName)
    {
      std::string::size_type extensionIndex =
        fileName.rfind(extensionDelimiter());
      if(extensionIndex == std::string::npos) {
        return std::make_pair(fileName, std::string(""));
      }
      std::string::size_type pathIndex =
        fileName.rfind(pathDelimiter());
      if((pathIndex != std::string::npos)
         && (pathIndex > extensionIndex)) {
        return std::make_pair(fileName, std::string(""));
      }
      return std::make_pair(fileName.substr(0, extensionIndex),
                            fileName.substr(extensionIndex, std::string::npos));
    }


    // This function accepts a path returns a pair of strings in which
    // the first element is the directory name and the second is the
    // filename.
    std::pair<std::string, std::string>
    splitPath(std::string const& path)
    {
      return brick::portability::splitPath(path);
    }

  } // namespace utilities

} // namespace brick
