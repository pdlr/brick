/**
***************************************************************************
* @file brick/common/Exception.hh
*
* Header file declaring some exception types.
*
* Copyright (c) 2003-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_EXCEPTION_HH
#define BRICK_COMMON_EXCEPTION_HH

#include <cstring>
#include <exception>


/**
 ** This macro to defines the maximum length of the exception internal
 ** error message string.  That is, it defines how long the what()
 ** message can be.
 **/
#ifndef BRICK_EXCEPTION_MESSAGE_LENGTH
#define BRICK_EXCEPTION_MESSAGE_LENGTH 512
#endif /* #ifndef BRICK_EXCEPTION_MESSAGE_LENGTH */


/**
 ** This macro controls how much extra storage will be statically
 ** allocated in each exception instance to handle additional
 ** user-defined payload data.  It specifies how much extra stuff you
 ** can piggyback onto a brick::exception if you want to communicate
 ** with calling contexts in a way that doesn't involve the "what"
 ** message.
 **
 ** WARNING: If you change the value of this macro, then the size of
 ** brick::Exception will change, and you must recompile all of your
 ** code with the new value or you'll have nasty linkage errors.
 **/
#ifndef BRICK_EXCEPTION_PAYLOAD_SIZE
#define BRICK_EXCEPTION_PAYLOAD_SIZE 8192
#endif /* #ifndef BRICK_EXCEPTION_PAYLOAD_SIZE */


/**
 * The BRICK_THROW macro constructs an exception instance using the
 * specified error message & function name, and automatically provides
 * file name & line number arguments.  You might use it like this:
 *
 *   BRICK_THROW(brick::IndexException, "myFunction(int)",
 *             "Argument value is out-of-bounds.");
 *
 * @param exception This argument specifies the type of exception to
 * throw.
 *
 * @param functionName This argument specifies the name of the
 * function from which the exception is being thrown.
 *
 * @param message This argument specifies a C-style string describing
 * the error.
 */
#define BRICK_THROW(exception, functionName, message) {         \
    throw exception(message, functionName, __FILE__, __LINE__); \
  }


/**
 * The BRICK_THROW2 macro constructs an exception instance using the
 * specified error message, and automatically provides file name &
 * line number arguments.  It differs from BRICK_THROW in that the
 * function name is not specified.
 *
 *   BRICK_THROW(brick::IndexException, "Bad state in graphics hardware.");
 *
 * @param exception This argument specifies the type of exception to
 * throw.
 *
 * @param message This argument specifies a C-style string describing
 * the error.
 */
#define BRICK_THROW2(exception, message) {              \
    throw exception(message, __FILE__, __LINE__);       \
  }


/**
 * This macro makes it easy to declare new exception types.  The point
 * of the Exception class heirarchy is to encourage client code to
 * derive its own exception classes.  Unfortunately, there are several
 * constructors to declare, which makes deriving exception subclasses
 * a headache, so nobody does it.  The maintainability hit we take by
 * having a magic macro for declaring exception types is less than the
 * maintainability hit we take by always throwing generic exceptions,
 * so we bite the bullet and define a macro for subclassing
 * exceptions.  For example, if you needed some subclasss of
 * brick::IOException to communicate parsing failures in a file-parsing
 * class, you might use it like this:
 *
 * @code
 * BRICK_DECLARE_EXCEPTION_TYPE(ParseException, brick::common::IOException);
 * BRICK_DECLARE_EXCEPTION_TYPE(SyntaxException, ParseException);
 * BRICK_DECLARE_EXCEPTION_TYPE(MissingArgumentException, ParseException);
 * @endcode
 *
 * @param ExceptionName This argument is the name of the new exception
 * class.
 *
 * @param ParentName This argument is the name of the parent class.
 * It can be set to brick::Exception or to any class derived from
 * brick::Exception using the BRICK_DECLARE_EXCEPTION_TYPE macro.
 */
#define BRICK_DECLARE_EXCEPTION_TYPE(ExceptionName, ParentName)         \
  class ExceptionName : public ParentName {                             \
  public:                                                               \
  ExceptionName() throw()                                               \
  : ParentName("", #ExceptionName) {}                                   \
                                                                        \
  ExceptionName(const char* message) throw()                            \
    : ParentName(message, #ExceptionName) {}                            \
                                                                        \
  ExceptionName(const char* message, const char* fileName,              \
                int lineNumber) throw()                                 \
    : ParentName(message, #ExceptionName, 0, fileName, lineNumber) {}   \
                                                                        \
  ExceptionName(const char* message, const char* functionName,          \
                const char* fileName, int lineNumber) throw()           \
    : ParentName(message, #ExceptionName, functionName, fileName,       \
                 lineNumber) {}                                         \
                                                                        \
  ExceptionName(const ExceptionName& source) throw()                    \
    : ParentName(source) {}                                             \
                                                                        \
  virtual ~ExceptionName() throw() {}                                   \
                                                                        \
  protected:                                                            \
  ExceptionName(const char* message, const char* childClassName) throw() \
    : ParentName(message, childClassName) {}                            \
                                                                        \
  ExceptionName(const char* message, const char* childClassName,        \
                const char* functionName, const char* fileName,         \
                int lineNumber) throw()                                 \
    : ParentName(message, childClassName, functionName, fileName,       \
                 lineNumber) {}                                         \
  }


/**
 ** This namespace comprises all of the symbols defined in the brick utility
 ** libraries.
 **/
namespace brick {

  /**
   ** This namespace contains classes, functions and typedefs for
   ** working with exceptions, reference counting, portable numeric
   ** types, checking and switching byte order, extending the standard
   ** library, and more.
   **/
  namespace common {

    /**
     ** Base class for all exceptions thrown from code in namespace brick.
     ** Note that if you must handle every possible exception
     ** from a routine in namespace crl, you need to catch the
     ** standard exceptions explicitly, like this:
     **
     ** @code
     **   try {
     **     brick::foo();
     **   } catch(brick::Exception& e) {
     **     // ...
     **   } catch {std::exception& e) {
     **     // ...
     **   }
     ** @endcode
     **
     ** if you must handle every possible exception from a routine in
     ** namespace brick, and you don't care about the added features of
     ** brick::Exception, you can omit the catch(brick::Exception&) block
     ** and still catch everything, like this:
     **
     ** @code
     **   try {
     **     brick::foo();
     **   } catch {std::exception& e) {
     **     // ...
     **   }
     ** @endcode
     **
     **
     ** Finally, note that brick::Exception and its immediate subclasses
     ** use no dynamically allocated memory.  Presumably this will pay
     ** off eventually.
     **/
    class Exception : public std::exception {
    public:

      /**
       * This constructor sets the internal "what()" message.  The
       * message should generally provide information about the test
       * failure.  For example: "Unable to acquire global lock file:
       * /var/myApp/.lockFile"
       *
       * @param message This argument is a C-style string specifying the
       * text of the message.
       */
      Exception(const char* message) throw();


      /**
       * This constructor builds the internal "what()" message using
       * detailed information about the error.
       *
       * @param message This argument specifies a description of the
       * failure.
       *
       * @param fileName This argument specifies the name of the file in
       * which the failure occurred.
       *
       * @param lineNumber This argument specifies the line number at
       * which the failure occurred.
       */
      Exception(const char* message, const char* fileName,
                int lineNumber) throw();


      /**
       * This constructor builds the internal "what()" message using
       * detailed information about the failure.
       *
       * @param message This argument specifies a description of the
       * failure.
       *
       * @param functionName This argument specifies the name of the
       * function in which the failure occurred.
       *
       * @param fileName This argument specifies the name of the file in
       * which the failure occurred.
       *
       * @param lineNumber This argument specifies the line number at
       * which the failure occurred.
       */
      Exception(const char* message, const char* functionName,
                const char* fileName, int lineNumber) throw();


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       */
      Exception(const Exception& source) throw();


      /**
       * Destructor.
       */
      virtual ~Exception() throw() {}


      /**
       * This public method copies user-supplied data out of the
       * exception class.  It is used in conjunction with setPayload()
       * to allow the throwing context to communicate with the
       * catching context without affecting the what message.
       *
       * @param buffer This argument must point to at least
       * BRICK_EXCEPTION_PAYLOAD_SIZE bytes of valid memory space.
       * Some of this space will be overwritten with a copy of the
       * exception's payload data.
       *
       * @param payloadSize This argument will be set to indicate how
       * many bytes of data were copied.
       */
      virtual void
      getPayload(char* buffer, unsigned int& payloadSize) const throw();


      /**
       * This public method copies returns the current length of the
       * user-supplied data stored in the exception class.  It is
       * useful if you need to append to that data using the three
       * argument version of member function setPayload().  You might
       * use it like this:
       *
       * @code
       *   char[] message = "To be appended.\n";
       *   myException.setPayload(myException.getPayloadSize(),
       *                          message, sizeof(message));
       * @endcode
       *
       * @param return The return value indicates the current size of
       * the user supplied payload data.
       */
      virtual unsigned int
      getPayloadSize() const throw() {return m_payloadSize;}


      /**
       * This public method copies user-supplied data into the
       * exception class.  It is used in conjunction with getPayload()
       * to allow the throwing context to communicate with the
       * catching context without affecting the what message.
       *
       * @param buffer This argument points to a data buffer from
       * which at most BRICK_EXCEPTION_PAYLOAD_SIZE bytes will be
       * copied.
       *
       * @param bufferSize This argument controls how many bytes are
       * copied from buffer.
       */
      virtual void
      setPayload(char const* buffer, unsigned int bufferSize) throw();


      /**
       * This public method copies user-supplied data into the
       * exception class.  It differs from the two-argument version of
       * setPayload() in that the user-supplied data is copied after
       * skipping over a user-determined number of bytes of internal
       * storage.  See member function getPayloadSize() for an example
       * of how to use this member function.  If the existing payload
       * is shorter than skipBytes, then any uninitialized bytes of
       * the payload will be filled in with \0.
       *
       * @param skipBytes This argument specifies how many bytes of
       * exisiting payload to skip before beginning to copy the input
       * buffer.
       *
       * @param buffer This argument points to a data buffer from
       * which at most BRICK_EXCEPTION_PAYLOAD_SIZE bytes will be
       * copied.
       *
       * @param bufferSize This argument controls how many bytes are
       * copied from buffer.
       */
      virtual void
      setPayload(unsigned int skipBytes, char const* buffer,
                 unsigned int bufferSize) throw();


      /**
       * This public method returns a C-style string describing the
       * condition which caused the exception to be thrown.  The value
       * of this string is set during construction based on the provided
       * constructor arguments.
       *
       * @return The return value is a C-style string describing the
       * condition which caused the exception to be thrown.
       */
      virtual const char* what() const throw() {return m_message;}


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument specifies the Exception instance to
       * be copied.
       *
       * @return The return value is a reference to *this.
       */
      virtual Exception& operator=(const Exception& source) throw();

    protected:
      /**
       * This protected constructor is provided so that derived classes
       * can easily personalize this->what() output. The what() message
       * will be constructed by combining the two arguments, separated
       * by the two character string ": ".  If the childClassName
       * argument is set to 0, the what() message will simply be copied
       * from message.  For example, calling the constructor like this:
       *
       * @code
       *   brick::Exception("My message.", "MyException")
       * @endcode
       *
       * results in the following what() message:
       *
       *   "MyException: My message."
       *
       * Calling the constructor like this:
       *
       * @code
       *   brick::Exception("My message.", 0)
       * @endcode
       *
       * results in the following what() message:
       *
       *   "My message."
       *
       * @param message This argument is a C-style string specifying the
       * desired what() message.
       *
       * @param childClassName This argument, if not equal to 0,
       * specifies the name of the child class from which this
       * constructor was called.  This name will be prepended to the
       * what message.
       */
      Exception(const char* message, const char* childClassName) throw();


      /**
       * This protected constructor is provided so that derived classes
       * can easily generate standardized this->what() output. The
       * what() message will be constructed by combining all of the
       * arguments in a visually pleasing way.  The constructor might be
       * called from a derived class like this:
       *
       * @code
       *   brick::Exception("My message.", "MyException", "myFunction(int)",
       *                    "myFile.cc", 265)
       * @endcode
       *
       * @param message This argument is a C-style string specifying the
       * desired what() message.
       *
       * @param childClassName This argument, if not equal to 0,
       * specifies the name of the child class from which this
       * constructor was called.  This name will be prepended to the
       * what message.
       *
       * @param functionName This argument specifies the name of the
       * function from which the exception was thrown.
       *
       * @param fileName This argument specifies the name of the source
       * file defining the function from which the exception was thrown.
       *
       * @param lineNumber This argument specifies the source line
       * number at which the exception was thrown.
       *
       */
      Exception(const char* message, const char* childClassName,
                const char* functionName, const char* fileName,
                int lineNumber) throw();


      /**
       * This protected member is a C-style string which holds a
       * message describing the condition that caused the exception to
       * be thrown.
       */
      char m_message[BRICK_EXCEPTION_MESSAGE_LENGTH];

      char m_payload[BRICK_EXCEPTION_PAYLOAD_SIZE];
      unsigned int m_payloadSize;
    };


    /**
     ** This is an Exception class for errors in which the outside world
     ** is misbehaving.  This is also an example of how to subclass
     ** Exception without using the BRICK_DECLARE_EXCEPTION_TYPE macro.
     **/
    class IOException : public Exception {
    public:
      IOException() throw()
        : Exception("", "IOException") {}

      IOException(const char* message) throw()
        : Exception(message, "IOException") {}

      IOException(const char* message, const char* fileName,
                  int lineNumber) throw()
        : Exception(message, "IOException", 0, fileName, lineNumber) {}

      IOException(const char* message, const char* functionName,
                  const char* fileName, int lineNumber) throw()
        : Exception(message, "IOException", functionName, fileName,
                    lineNumber) {}

      IOException(const IOException& source) throw()
        : Exception(source) {}

      virtual ~IOException() throw() {}

    protected:
      IOException(const char* message, const char* childClassName) throw()
        : Exception(message, childClassName) {}

      IOException(const char* message, const char* childClassName,
                  const char* functionName, const char* fileName,
                  int lineNumber) throw()
        : Exception(message, childClassName, functionName, fileName,
                    lineNumber) {}
    };


    /**
     ** This is an Exception class for errors in which an array index or
     ** similar argument is out of bounds.
     **/
    class IndexException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(IndexException, Exception);


    /**
     ** This is an Exception class for errors which could have been
     ** caught at compile time.
     **/
    class LogicException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(LogicException, Exception);


    /**
     ** This type of Exception is thrown when a piece of code has not
     ** been written, and the developer had the foresight to document
     ** its non-existence.
     **/
    class NotImplementedException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(NotImplementedException, Exception);


    /**
     ** This exception is thrown when an error occurs which could not
     ** have been anticipated at compile time, and for which ValueException,
     ** IOException, etc., are not appropriate.
     **/
    class RunTimeException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(RunTimeException, Exception);


    /**
     ** This exception is thrown when the internal state of a class is
     ** inconsistent, or when the calling environment interacts with a class
     ** in a way which is inconsistent with its internal state.
     **/
    class StateException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(StateException, Exception);


    /**
     ** This exception is thrown when the argument to a function has
     ** and inappropriate value, and when a more specific exception
     ** is not appropriate.
     **/
    class ValueException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(ValueException, Exception);

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_EXCEPTION_HHH */
