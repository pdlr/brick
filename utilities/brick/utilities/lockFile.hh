/**
***************************************************************************
* @file brick/utilities/lockFile.hh
*
* Header file declaring the LockFile class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_LOCKFILE_HH
#define BRICK_UTILITIES_LOCKFILE_HH

#include <string>

namespace brick {

  namespace utilities {
    
    /**
     ** The LockFile class tries to create a file that is uniquely
     ** owned by the calling process.  If a file with the specified
     ** name already exists, then it will not be created, and
     ** subsequent calls to the isValid() method of the LockFile
     ** instance will return false.  If the file does not exist, and
     ** is succesfully created, then subsequent calls to the isValid()
     ** method will return true, and the file will be deleted when the
     ** LockFile instance is destroyed.
     **
     ** NOTE: This class currently does not compile under windows.
     **
     ** WARNING: The file creation method used by LockFile is not
     ** atomic on NFS filesystems.  This means that using LockFile on
     ** an NFS filesystem is likely to introduce a race condition to
     ** your software.
     **
     ** If you want your code to block waiting for a lock, you might
     ** use LockFile like this:
     **
     ** @code
     ** {
     **   // This constructor call will block until the lock is obtained.
     **   LockFile lockFile("/var/tmp/myLockDir/0000.lock", "message", -1);
     **   myFunctionThatRequiresLocking();
     ** } // LockFile destructor releases lock here.
     ** @endcode
     ** 
     ** Here's another way to use LockFile:
     **
     ** @code
     ** {
     **   // This constructor call will return immediately.
     **   LockFile lockFile("/var/tmp/myLockDir/0000.lock", "message");
     **   while(!lockFile.retry(1.0)) {
     **     std::cout << "Still waiting for lock..." << std::endl;
     **   }
     **   myFunctionThatRequiresLocking();
     ** } // LockFile destructor releases lock here.
     ** @endcode
     ** 
     **/
    class LockFile {
    public:

      /** 
       * The constructor attempts to create the lock file.
       * 
       * @param fileName This argument specifies the file to be created.
       *
       * @param timeout The constructor will wait this long (in
       * seconds) trying to establish the lock.  If the attempt to
       * acquire the lock is still unsuccessful by the end of the
       * timeout (perhaps because another LockFile instance has
       * already created a lock on this filename) the constructor will
       * complete execution, but this->isValid() will return false.
       * Setting this argument to 0.0 will cause the constructor to
       * return immediately.  Setting this argument less than zero
       * allow the constructor to block indefinitely waiting for the
       * lock.
       */
      LockFile(std::string const& fileName,
               double timeout = 0.0);


      /** 
       * The constructor attempts to create a non-empty lock file.  If
       * the lock file is successfully created, then the specified
       * contents string will be written to the file.
       * 
       * @param fileName This argument specifies the file to be created.
       *
       * @param contents This string is to be written to the lock file
       * after creation.  You might use it to associate a useful message
       * with the file, such as "This lock file created by process
       * number 2412."
       *
       * @param timeout The constructor will wait this long (in
       * seconds) trying to establish the lock.  If the attempt to
       * acquire the lock is still unsuccessful by the end of the
       * timeout (perhaps because another LockFile instance has
       * already created a lock on this filename) the constructor will
       * complete execution, but this->isValid() will return false.
       * Setting this argument to 0.0 will cause the constructor to
       * return immediately.  Setting this argument less than zero
       * allow the constructor to block indefinitely waiting for the
       * lock.
       */
      LockFile(std::string const& fileName, std::string const& contents,
               double timeout = 0.0);


      /**
       * The destructor destroys the LockFile instance and deletes any
       * file created by the constructor, thereby releasing the lock.
       */
      ~LockFile();


      /** 
       * This member function reports whether or not the lock was
       * successfully obtained.  That is, it reports whether or not
       * the constructor or retry() method was successful in asserting
       * the lock.
       * 
       * @return The return value is true if the lock was obtained,
       * false otherwise.
       */
      bool
      isValid();


      /** 
       * This member function tries to obtain the lock, and is useful
       * for cases where the constructor failed to obtain the lock.
       * Calling this function when the file has already been locked
       * has no effect.
       * 
       * @param timeout This argument controls how long the call to
       * retry() is allowed to wait for a lock to be established.  Its
       * meaning is just like that of the identically named
       * constructor argument.
       *
       * @return The return value is true if the lock was obtained,
       * false otherwise.
       */
      bool
      retry(double timeout = 0.0);
      
    private:

      int m_fileDescriptor;
      std::string m_fileName;

      // This member is currently only to soothe compiler warnings.
      int m_returnCode;
      
    }; // class LockFile

  } // namespace utilities
  
} // namespace brick

#endif /* #ifndef BRICK_UTILITES_LOCKFILE_HH */
