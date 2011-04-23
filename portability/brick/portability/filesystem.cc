/**
***************************************************************************
* @file brick/portability/filesystem.cpp
* 
* Source file defining portability routines for dealing with filesystems.
*
* Copyright (C) 2003-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/portability/filesystem.hh>

/* ===================== Common includes ===================== */

#include <algorithm>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <brick/common/exception.hh>

/* ===================== End common includes ===================== */

#ifdef _WIN32

/* ===================== Windows includes ===================== */

#include <io.h>
#include <windows.h>

/* ===================== End Windows includes ===================== */

#else /* #ifdef _WIN32 */

/* ===================== Linux includes ===================== */

#include <dirent.h>
#include <errno.h>
#include <unistd.h>

/* ===================== End Linux includes ===================== */

#endif /* #ifdef _WIN32 */


/* ===================== Common code ===================== */

namespace brick {

  namespace portability {


    // This function returns a bool indicating whether or not the
    // specified path is a directory.
    bool
    isDirectory(const std::string& path)
    {
      struct stat statBuf;
      memset(&statBuf, 0, sizeof(statBuf));
      int returnValue = stat(path.c_str(), &statBuf);
      if(returnValue == 0) {
        // if(S_ISDIR(statBuf.st_mode)) {
        if((statBuf.st_mode & S_IFMT) == S_IFDIR) {
          return true;
        }
      }
      return false;
    }
    

    // Joins two path elements with the appropriate delimiter.
    std::string
    joinPath(const std::string& part0, const std::string& part1)
    {
      if(part0.empty()) {
        return part1;
      }
      if(part1.empty()) {
        return part0;
      }
      // If part1 is an absolute path, ignore part0.
      if(part1.size() >= pathDelimiter().size()) {
        if(part1.compare(0, pathDelimiter().size(), pathDelimiter()) == 0) {
          return part1;
        }
      }
      // If part0 ends in a path delimiter, simply concatenate the two parts.
      if(part0.size() >= pathDelimiter().size()) {
        if(part0.compare(
             part0.size() - pathDelimiter().size(), std::string::npos,
             pathDelimiter()) == 0) {
          return part0 + part1;
        }
      }
      // Add the missing delimiter and return.
      return part0 + pathDelimiter() + part1;
    }


    // This function accepts a path returns a pair of strings in which
    // the first element is the directory name and the second is the
    // filename.
    std::pair<std::string, std::string>
    splitPath(const std::string& path)
    {
      typedef std::vector<std::string>::const_iterator DelimIter;
      
      const std::vector<std::string>& delimiterVector = pathDelimiters();
      std::string::size_type delimiterIndex = std::string::npos;

      DelimIter delimIter = delimiterVector.begin();
      while(delimIter != delimiterVector.end()) {
        delimiterIndex = std::min(delimiterIndex, path.rfind(*delimIter));
        ++delimIter;
      }
      if(delimiterIndex == std::string::npos) {
        return std::make_pair(std::string(""), path);
      }
      if(delimiterIndex == path.size() - 1) {
        return std::make_pair(path, std::string(""));
      }
      return std::make_pair(
        path.substr(0, delimiterIndex + 1),
        path.substr(delimiterIndex + 1, std::string::npos));
    }
    
  } // namespace portability

} // namespace brick

/* ===================== End common code ===================== */


#ifdef _WIN32

/* ===================== Windows code ===================== */

// Anonymous namespace for stuff local to this file.
namespace brick {

  namespace portability {

    // Initialize on first use to avoid static initialization order fiasco.
    const std::string& pathDelimiter() {
      static std::string delimiter = "\\";
      return delimiter;
    }


    // Initialize on first use to avoid static initialization order fiasco.
    const std::vector<std::string>& pathDelimiters() {
      static bool needsInit = true;
      static std::vector<std::string> delimiterVector;
      if(needsInit) {
        delimiterVector.push_back("\\");
        delimiterVector.push_back("/");
        needsInit = false;
      }
      return delimiterVector;
    }

    
    // Initialize on first use to avoid static initialization order fiasco.
    const std::string& extensionDelimiter() {
      static std::string delimiter = ".";
      return delimiter;
    }


    // Returns the names of the entries in the specified directory, in
    // no particular order.
    std::vector<std::string>
    listDirectory(const std::string& directoryName, bool fullPath)
    {
      std::vector<std::string> listing;
      std::string fileNameGlob = joinPath(directoryName, "*");
      WIN32_FIND_DATA findData;
      HANDLE handle = FindFirstFile(fileNameGlob.c_str(), &findData); 
      if (handle == INVALID_HHANDLE_VALUE) { 
        return listing;
      } 
      while(1) {
        listing.push_back(std::string(findData.cFileName));
        if (!FindNextFile(handle, &findData)) {
          if (GetLastError() == ERROR_NO_MORE_FILES) {
            break;
          }
          std::ostringstream message;
          message << "Problem reading file names from " << directoryName;
          BRICK_THROW(brick::common::IOException, "listDirectory()",
                      message.str().c_str());
        }
      }
      if (!FindClose(handle)) { 
        std::ostringstream message;
        message << "Problem closing search handle for " << directoryName;
        BRICK_THROW(brick::common::IOException, "listDirectory()",
                    message.str().c_str());
      }
    
      if(fullPath) {
        // std::transform(listing.begin(), listing.end(), listing.begin(),
        //                std::bind1st(std::ptr_fun(joinPath), directoryName));
        typedef std::vector<std::string>::iterator Iter;
        for(Iter iter = listing.begin(); iter != listing.end(); ++iter) {
          *iter = joinPath(directoryName, *iter);
        }
      }
      return listing;
    }
  
  } // namespace portability
  
} // namespace brick

/* ===================== End Windows code ===================== */

#else /* #ifdef _WIN32 */

/* ===================== Linux code ===================== */

namespace brick {

  namespace portability {

    // Initialize on first use to avoid static initialization order fiasco.
    const std::string& pathDelimiter() {
      static std::string delimiter = "/";
      return delimiter;
    }

  
    // Initialize on first use to avoid static initialization order fiasco.
    const std::vector<std::string>& pathDelimiters() {
      bool needsInit = true;
      static std::vector<std::string> delimiterVector;
      if(needsInit) {
        delimiterVector.push_back("/");
        needsInit = false;
      }
      return delimiterVector;
    }

    
    // Initialize on first use to avoid static initialization order fiasco.
    const std::string& extensionDelimiter() {
      static std::string delimiter = ".";
      return delimiter;
    }


    // Returns the names of the entries in the specified directory, in
    // no particular order.
    std::vector<std::string>
    listDirectory(const std::string& directoryName, bool fullPath)
    {
      std::vector<std::string> listing;

      DIR* directoryPtr = opendir(directoryName.c_str());
      if(directoryPtr == 0) {
        std::ostringstream message;
        char charBuffer[1024];
        if(strerror_r(errno, charBuffer, 1024) == 0) {
          charBuffer[1023] = '\0';
          message << charBuffer;
        } else {
          message << "Unknown error opening directory";
        }
        message << ": " << directoryName << std::endl;
        BRICK_THROW(brick::common::IOException, "listDirectory()",
                    message.str().c_str());
      }
      try {
        while(1) {
          struct dirent* directoryEntryPtr = readdir(directoryPtr);
          if(directoryEntryPtr == 0) {
            // NULL return value indicates no more entries.
            // Can also indicate a bat directory stream, but we just
            // checked that above.
            break;
          }
          listing.push_back(std::string(directoryEntryPtr->d_name));
        }
      } catch(...) {
        closedir(directoryPtr);
        throw;
      }
      closedir(directoryPtr);
      
      if(fullPath) {
        // std::transform(listing.begin(), listing.end(), listing.begin(),
        //                std::bind1st(std::ptr_fun(joinPath), directoryName));
        typedef std::vector<std::string>::iterator Iter;
        for(Iter iter = listing.begin(); iter != listing.end(); ++iter) {
          *iter = joinPath(directoryName, *iter);
        }
      }
      return listing;
    }

  } // namespace portability

} // namespace brick

/* ===================== End Linux code ===================== */

#endif /* #ifdef _WIN32 */

