/**
***************************************************************************
* @file brick/utilities/lockFile.cc
*
* Source file defining the LockFile class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <brick/portability/standardC.hh>
#include <brick/utilities/lockFile.hh>
#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace utilities {
    
    // The constructor attempts to create the lock file.
    LockFile::
    LockFile(std::string const& fileName,
             double timeout)
      : m_fileDescriptor(open(fileName.c_str(), O_WRONLY | O_CREAT | O_EXCL,
                              S_IRUSR)),
        m_fileName(fileName)
    {
      if(timeout != 0.0) {
        double stopTime = std::numeric_limits<double>::max();
        double sleepTime = 0.1;
        if(timeout > 0.0) {
          stopTime = getCurrentTime() + timeout;
          sleepTime = std::min(1.0, timeout / 100.0);
        }
        while(m_fileDescriptor < 0) {
          portableSleep(sleepTime);
          m_fileDescriptor =
            open(fileName.c_str(), O_WRONLY | O_CREAT | O_EXCL, S_IRUSR);
          if(getCurrentTime() > stopTime) {
            break;
          }
        }
      }
      if(m_fileDescriptor != -1) {
        close(m_fileDescriptor);
      }
    }
  

    // The constructor attempts to create a non-empty lock file.
    LockFile::
    LockFile(std::string const& fileName,
             std::string const& contents,
             double timeout)
      : m_fileDescriptor(open(fileName.c_str(), O_WRONLY | O_CREAT | O_EXCL,
                              S_IWUSR)),
        m_fileName(fileName)
    {
      if(timeout != 0.0) {
        double stopTime = std::numeric_limits<double>::max();
        double sleepTime = 0.1;
        if(timeout > 0.0) {
          stopTime = getCurrentTime() + timeout;
          sleepTime = std::min(1.0, timeout / 100.0);
        }
        while(m_fileDescriptor < 0) {
          portableSleep(sleepTime);
          m_fileDescriptor =
            open(fileName.c_str(), O_WRONLY | O_CREAT | O_EXCL, S_IRUSR);
          if(getCurrentTime() > stopTime) {
            break;
          }
        }
      }
      if(m_fileDescriptor != -1) {
        // Success of these function calls is largely irrelevant.
        write(m_fileDescriptor, (const void*)contents.c_str(), contents.size());
        close(m_fileDescriptor);
        chmod(m_fileName.c_str(), S_IRUSR);
      }
    }
      

    // The destructor destroys the LockFile instance and deletes any
    // file created by the constructor, thereby releasing the lock.
    LockFile::
    ~LockFile()
    {
      if(m_fileDescriptor != -1) {
        // There's no really good way to handle a failure in one of
        // these two calls, since we don't want to throw an exception
        // from a destructor.  Instead we just assume they'll succeed.
        chmod(m_fileName.c_str(), S_IWUSR);
        unlink(m_fileName.c_str());
      }
    }


    // This method reports whether or not the lock was successfully
    // obtained.
    bool
    LockFile::
    isValid()
    {
      return (m_fileDescriptor != -1);
    }

  } // namespace utilities
  
} // namespace brick
